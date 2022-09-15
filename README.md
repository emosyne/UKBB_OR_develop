# 2HH
# Emanuele's first workflow
eosimo@ic.ac.uk

Run this example pipeline with:

```
nextflow run main.nf -profile servername
```

To install additional modules, use the [nf-core tools](https://github.com/nf-core/tools), for example:

```
nf-core modules install star/align
```

To create your own module, also use the [nf-core tools](https://github.com/nf-core/tools):

```
nf-core modules create
```

It will ask you for the module name and other parameters.

After creating or installing a module, you need to add an `include` instruction for that module into your workflow and also add a corresponding structure to `modules.config` (see `workflows/example_wf.nf` and `modules.config` for examples).
