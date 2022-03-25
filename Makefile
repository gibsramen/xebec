TMPDIR := $(shell mktemp -d)
COOKIE_DIR := $(shell pwd)
TABLE_FILE := $(shell realpath tests/data/table.biom)
MD_FILE := $(shell realpath tests/data/metadata.tsv)
TREE_FILE := $(shell realpath tests/data/tree.tre)

all: pytest snaketest

test:
	pytest

snaketest:
	@cd $(TMPDIR); \
	cookiecutter --no-input $(COOKIE_DIR) feature_table_file=$(TABLE_FILE) sample_metadata_file=$(MD_FILE) phylogenetic_tree_file=$(TREE_FILE); \
	cd diversity-benchmark; \
	snakemake --cores 1; \
	rm -r $(TMPDIR)
