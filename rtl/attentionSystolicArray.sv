`timescale 1ns / 1ps
`default_nettype none


module attentionSystolicArray #(
    parameter int unsigned N = 4,
    parameter int unsigned d = 4,
    parameter int unsigned K = 4
)(
    input logic clk,
    input logic reset,
    input logic valid_input,
    input real K_matrix [0:N-1][0:d-1],  // K matrix input
    input real factorial_arr [0:K],  // Weight input for each PE
    input real Q_matrix   [0:N-1][0:d-1],  // Q matrix input
    input real V_matrix  [0:N-1][0:d-1],  // V matrix input
    input real softmax_temp, // Softmax temp (1/sqrt(d))
    /* verilator lint_off UNDRIVEN */
    output real Q_mult_K[0:N-1][0:N-1],  // Result output
    /* verilator lint_off UNDRIVEN */
    output real exp_Q_mult_K[0:N-1][0:N-1], // Result exponentiation of matrix
    /* verilator lint_off UNDRIVEN */
    output real exp_K_mult_Q_mult_V[0:N-1][0:d-1],  // Result output
    output real attention[0:N-1][0:d-1], // Result output
    output logic valid_result
);

    // Control Counter Logic
    localparam int unsigned MULT_CYCLES = 3*N+1+K+d;

    int unsigned counter_q, counter_d;

    always_ff @(posedge clk) begin
        if (reset)
            counter_q <= 0;
        else
            counter_q <= counter_d;
    end

    always_comb begin
        if (doProcess_q)
            counter_d = counter_q + 1;
        else
            counter_d = 0;
    end

    // Valid Result Signal Logic
    logic valid_result_q;

    always_ff @(posedge clk) begin
        if (reset)
            valid_result_q <= 0;
        else if (counter_q == MULT_CYCLES)
            valid_result_q <= 1;
        else
            valid_result_q <= 0;
    end

    assign valid_result = valid_result_q;

    // Process Control Logic
    logic doProcess_d, doProcess_q;

    always_ff @(posedge clk) begin
        if (reset)
            doProcess_q <= 0;
        else
            doProcess_q <= doProcess_d;
    end

    always_comb begin
        if (valid_input && weight_load_done)
            doProcess_d = 1;
        else if (counter_q == MULT_CYCLES + 1)
            doProcess_d = 0;
        else
            doProcess_d = doProcess_q;
    end

    // Weight Load Control Logic
    logic weight_load_enable_q;
    logic weight_load_done;
    logic [1:0] weight_load_counter;

    always_ff @(posedge clk) begin
        if (reset) begin
            weight_load_enable_q <= 1;
            weight_load_counter <= 0;
        end else if (weight_load_enable_q) begin
            weight_load_counter <= weight_load_counter + 1;
            if (weight_load_counter == 2) // Adjust based on loading requirements
                weight_load_enable_q <= 0;
        end
    end

    assign weight_load_done = !weight_load_enable_q;

    // Data Input Array
    real data_in_array [0:d-1];

    // Generate data_in_array based on the current counter_q value
    always_ff @(posedge clk) begin
      if (reset) begin
        for (int j = 0; j < d; j++) begin
            data_in_array[j] <= '0;
        end
      end else if (doProcess_d) begin
        for (int j = 0; j < d; j++) begin
          if ($signed(counter_d) - $signed(j) >= 0 && $signed(counter_d) - $signed(j) < N) begin
              // Skew Q_matrix and read each row into data_in_array
              data_in_array[j] <= Q_matrix[counter_d - j][j];
          end else begin
              // Add 0 padding
              data_in_array[j] <= '0;
          end
        end
      end
    end

    real K_matrix_in [0:N-1][0:d-1];
    // K matrix is loaded into weights upside down
    always_comb begin
        for (int i = 0; i < N; i++) begin
            for (int j = 0; j < d; j++) begin
                K_matrix_in[i][j] = K_matrix[N-1-i][j];
            end
        end
    end

    // Intermediate result of QK^T
    real Q_mult_K_result [0:N-1];
    // Instantiate the weight-stationary systolic array
    systolic_array #(
        .N(N),
        .M(d)
    ) K_mult_Q_systolicArray (
        .clk(clk),
        .reset(reset),
        .weight_load_enable(weight_load_enable_q),
        .doProcess(doProcess_q),
        .weight_in(K_matrix_in),
        .data_in(data_in_array),
        .result_out(Q_mult_K_result)
    );

    // Intermediate result of exp(QK^T)
    real exp_Q_mult_K_result [0:N-1];
    // Row-wise normalizer in the softmax function
    real softmax_norm;

    real Q_mult_K_result_normed [0:N-1];
    always_comb begin
        for (int i = 0; i < N; i++) begin
            Q_mult_K_result_normed[i] = Q_mult_K_result[i] * softmax_temp;
        end
    end

    systolic_array_exp #(
        .K(K),  // Number of Taylor terms (columns)
        .N(N)   // Number of rows
    ) exp_systolic_array (
        .clk(clk),
        .reset(reset),
        .doProcess(doProcess_q),
        .factorial_arr(factorial_arr), // Factorial terms for each column
        .data_in(Q_mult_K_result_normed), // Input column vector
        .exp_out(exp_Q_mult_K_result),           // Output exponent approximation
        .exp_sum_out(softmax_norm)    // Output row-wise accumulated sum
    );

    // Corner turning before V matrix multiplication
    real data_in_v [0:N-1];
    always_comb begin
        for (int i = 0; i < N; i++) begin
            data_in_v[i] = exp_Q_mult_K_result[N-1-i];
        end
    end

    // Weight matrix for V matrix multiplication (must be rotated 90 degrees for correctness)
    real weight_V_matrix [0:d-1][0:N-1];
    always_comb begin
        for (int i = 0; i < d; i++) begin
            for (int j = 0; j < N; j++) begin
                weight_V_matrix[i][j] = V_matrix[j][d-1-i]; // Rotate 90 degrees counterclockwise
            end
        end
    end
    // Intermediate result of exp(QK^T)V
    real exp_Q_mult_K_mult_V_result [0:d-1];
    systolic_array #(
        .N(d),
        .M(N)
    ) V_mult_systolicArray (
        .clk(clk),
        .reset(reset),
        .weight_load_enable(weight_load_enable_q),
        .doProcess(doProcess_q),
        .weight_in(weight_V_matrix),
        .data_in(data_in_v),
        .result_out(exp_Q_mult_K_mult_V_result)
    );

    // Final result of Attention(Q, K, V)
    real attention_result [0:d-1];
    // Systolic array for row-wise normalization
    systolic_array_norm #(
        .N(d)
    ) softmax_norm_systolic_array (
        .clk(clk),
        .reset(reset),
        .doProcess(doProcess_q),
        .norm_in(softmax_norm),
        .data_in(exp_Q_mult_K_mult_V_result),
        .data_normed_out(attention_result)
    );

    // Collect intermediate result of QK^T
    always_ff @(posedge clk) begin
        for (int i = 0; i < N; i++) begin
            for (int j = 0; j < N; j++) begin
                if (counter_q == d + i + j) begin
                    Q_mult_K[j][i] <= Q_mult_K_result[N-1-i]; //exp_out[N-1-i];
                end
            end
        end
    end

    // Collect intermeditate result of exp(QK^T)
    always_ff @(posedge clk) begin
        for (int i = 0; i < N; i++) begin
            for (int j = 0; j < N; j++) begin
                if (counter_q == d + i + j + K) begin
                    exp_Q_mult_K[j][i] <= exp_Q_mult_K_result[N-1-i];
                end
            end
        end
    end
    
    // Collect intermediate results of exp(QK^T)V
    always_ff @(posedge clk) begin
        for (int i = 0; i < d; i++) begin
            for (int j = 0; j < N; j++) begin
                if (counter_q == d + N + i + j + K) begin
                    exp_K_mult_Q_mult_V[j][i] <= exp_Q_mult_K_mult_V_result[d-1-i];
                end
            end
        end
    end
    // Collect intermediate results of Attention(Q, K, V)
    always_ff @(posedge clk) begin
        for (int i = 0; i < d; i++) begin
            for (int j = 0; j < N; j++) begin
                if (counter_q == d + N + i + j + K + 1) begin
                    attention[j][i] <= attention_result[d-1-i];
                end
            end
        end
    end



endmodule

`resetall
