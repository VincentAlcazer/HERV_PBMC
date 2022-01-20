rule sample_complete:
    output:
        touch("results/completed/{s}_{method}_completed.txt")
    input:
        rules.merge_matrices.output

