#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

include { UKBB_OR_develop } from './workflows/UKBB_OR_develop.nf'

workflow {
    UKBB_OR_develop()
}