#!/usr/bin/env nextflow

/*
    Function to get the directory for a given entry
*/
def getEntryDir(entry) {
    // Function to get the directory for a given entry
    // AF-1000000000000001 -> 1000/0000/0000/0001
    return entry.replaceFirst(/^AF-/, '').replaceAll(/(\d{4})(?=\d)/, '$1/')
}
