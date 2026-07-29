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

## `annovar_intervar.sif`
ANNOVAR is registration-gated. You must obtain the tarball yourself from
<https://annovar.openbioinformatics.org/> and stage it before the build:

```bash
cd containers/
mkdir -p _staging
cp /path/to/annovar.latest.tar.gz _staging/annovar.tar.gz
singularity build annovar_intervar.sif annovar_intervar.def
```

The container ships tools only; the ANNOVAR `humandb` and InterVar `intervardb`
data directories are **bind-mounted at runtime** via `params.annovar_humandb`
and `params.intervar_dir`. Download them with
`setup/setup_annovar_intervar.sh` after the image is built.

## Do not push
Both images bundle licensed / semi-licensed content. Keep them on
`/sc/arion/work/$USER/singularity_cache/` — do not push to a public registry.
