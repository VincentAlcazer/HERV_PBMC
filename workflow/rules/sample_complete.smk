rule sample_complete:
    output:
        touch("results/completed/{s}_completed.txt")
    input:
        rules.stellarscope_pseudobulk.output,
        rules.stellarscope_individual.output

