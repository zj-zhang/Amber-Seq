PROJECT="test-run" snakemake --dag  | dot -Tsvg > img/test-run.svg
PROJECT="test-run" snakemake --dag  | dot -Tpng > img/test-run.png
