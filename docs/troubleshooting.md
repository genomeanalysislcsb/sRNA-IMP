# Troubleshooting

## Host indexes rebuild unexpectedly

Host Bowtie indexes are reused only when the expected `.ebwt` or `.ebwtl` files already exist next to the original host FASTA files.

## GitHub SSH push fails

Use the working SSH key explicitly when needed:

```bash
GIT_SSH_COMMAND='ssh -i ~/.ssh/aion -o IdentitiesOnly=yes' git push
```

## Resume does not skip a stage

A stage is skipped only if:
- the completion marker exists
- the expected output files for that stage exist

## Assembly or novel steps fail immediately

Check that both sample-specific paths are present:
- `assembly_fasta`
- `assembly_gff`
