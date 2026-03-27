#!/usr/bin/env bash

# Example paths for database download and setup
# Copy and edit as needed.

THREADS=32
READ_LENGTH=50

DB_ROOT="${PWD}/databases"

HUMAN_DIR="${DB_ROOT}/human"
NONHUMAN_DIR="${DB_ROOT}/nonhuman"
SHARED_DIR="${DB_ROOT}/shared"

RFAM_DIR="${SHARED_DIR}/rfam"

KRAKEN_DIR="${NONHUMAN_DIR}/kraken2/standard"
BRACKEN_DIR="${NONHUMAN_DIR}/bracken"

HUMAN_MIRNA_DIR="${HUMAN_DIR}/mirna"
HUMAN_TRNA_DIR="${HUMAN_DIR}/trna"
HUMAN_RRNA_DIR="${HUMAN_DIR}/rrna"
HUMAN_OTHERRNA_DIR="${HUMAN_DIR}/otherrna"

NONHUMAN_MIRNA_DIR="${NONHUMAN_DIR}/mirna"
NONHUMAN_TRNA_DIR="${NONHUMAN_DIR}/trna"
