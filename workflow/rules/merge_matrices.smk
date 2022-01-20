rule merge_matrices:
    conda:
        "../envs/merger.yaml"
    output:
        "results/counts_matrix_R/{s}/{s}_{method}_matrix.Rds"
    input:
        protein_coding_dir = "results/starsolo_alignment/{s}/{s}.Solo.out/Gene/filtered/",
        transposable_elements_dir = "results/telescope_{method}/{s}/"
    benchmark:
        "benchmarks/merge_matrices/{s}_{method}_merge_matrices.tsv"
    log:
        "results/counts_matrix_R/{s}/{s}_{method}.log"
    params:
        sample = "{s}_{method}"
    script:
        "../scripts/merge_count_matrices.R"

