### Primer File Information
These are the fasta files called in Cutadapt to search for primers to reomove. For 5' primers, the sequence is preceded by an ^. This ^ indicates that the primer is anchored, which means it is found at the 5' end of the sequence and is required to be found in full (i.e. the read will not be kept if the primer is not found).  Most of these primers (except RC primers, see below) have spacers added to the 5' end to increase basepair heterogeneity. These spacers need to be included in the primer sequence to be found, since these primers are anchored. These spacers follow those used and supplied by the Travis Glenn lab.

For primer files that end in "RC.fas", these are the reverse complement of their respective primers. These primers are not anchored (therefore the read is kept regardless of whether the primer is found). Cutadapt will search for these primers on the 3' end of the complementary read (i.e. the RC Forward primer will be found on the 3' end of the R2 (reverse) read). These are used to remove primers when there is read-through (i.e. when the amplicon is shorter than sequencing length).


We currently have more than eight primers pairs avalailable. To see a list of what is available, from the terminal (not the R Consule) use this script:
```
basename -a primers/*.fas | grep -v '_RC' | sed -E 's/-(F|R)\.fas$//' | sort -u
```

### Custom Primers
To make your own primer pairs for cutadapt to remove, follow the format of existing primers. The primer pair files should be named simply, ending in "-F.fas" and "-R.fas". Try to be as precise as possible in your primer name, so there is no ambiguity in what primers are being used. In the file, the primers need to have the same name as the primer filename preceeding "F.fas". If you are using primers that have spacers, your file needs to include the separate primer sequences for each spacer (and include the primer without spacers). All primer sequences need to have the same name. For example, MiFish-F.fas has five primer sequences, and all are named "MiFish". All resulting products from this pipeline will contain the name of the primer in the file. The primer names in the -F and -R files should be the same. Each sequence in the file should be preceded by a "^", because this is a 5' sequence that is required to exist and be at the start of the read.

If you are using a short gene region that may have read-through (i.e. the sequencer will read the primer/adapter on the 3' end), you should also include a reverse-complement primer file. This primer is not required, and may not be at the end of the read, so no need for a "^" before the sequence. You also don't need to include spacers. Only the actual primer sequence is needed. However, the sequence needs to be the reverse complement of the normal primer sequence. These files names should be the same as the 5' primer files, but end with "-F_RC.fas" and "-R_RC.fas". The name of the primers inside these files is not as important as in the 5' primer files.

Custom primers should be saved in the primer folder, along with the rest of the primers.

