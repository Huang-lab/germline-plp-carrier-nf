# Container images

These `.def` files are built **on Minerva**. Do not build them here (the
authoring environment); do not push either image to a public registry.

## `annotate.sif`
Freely-licensed tooling only.

```bash
ml singularity-ce
cd containers/
singularity build annotate.sif annotate.def
```

## ACMG (fastVEP) — no container
ACMG-AMP is produced by **fastVEP** (Huang-lab/fastVEP), a native Rust binary
built once via `acmg_fastvep/setup_fastvep.sh` — there is no ACMG container. The
`ACMG_FASTVEP` process runs the `fastvep` binary directly (on PATH), or in a
container if you set `params.container_fastvep`. See `acmg_fastvep/README.md`.

## Do not push
`annotate.sif` bundles freely-licensed tooling but keep it on
`/sc/arion/work/$USER/singularity_cache/` — do not push to a public registry.
