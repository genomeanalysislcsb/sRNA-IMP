# Workflow

![sRNA-IMP workflow](assets/workflow.svg)

## Resume Logic

Each workflow writes a completion marker under the sample directory. Re-runs skip completed stages when the marker and expected outputs are both present.
