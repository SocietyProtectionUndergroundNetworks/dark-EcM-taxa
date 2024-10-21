# dark-taxa
Data accompanying dark taxa paper

## Estimating the number of named EcM fungal species in the Catalog of Life

### Software needed:

1. [csvtk](https://bioinf.shenwei.me/csvtk/)
    * Follow the installation steps
    * Create an alias `tsvtk` pointing to the `csvtk` binary as that runs `csvtk` with the defaults for parsing and generating tab-separated-value (TSV) files rather than comma-separated-value (CSV) files
2. Unix command line with GNU Coreutils and perl installed
    * The code here was run on Ubuntu linux with perl pre-installed for generating a regexp, but any other string manipulation tool can also be used
    * Similarly - `cat`, `cut`, `sort`, and `uniq` commands were used from the GNU Coreutils package, but these steps can be done using any other standard variants on other OSs as well.

### Files needed:

1. [ECM_genera2.csv](ECM_genera2.csv) - List of 327 Ectomycorrhizal genera - one word/name per line
2. [NameUsage.tsv](NameUsage.tsv) - List of all Fungi names (accepted and synonym) downloaded from [Catalogue of Life](https://checklistbank.org)

Both files are included in this repo.

### Steps

1. Update NameUsage.tsv (Optional). If you want the latest version of this file, download all accepted fungal species names and synonyms from Catalogue of Life:

    * Go to [Catalogue of Life, Release 2024-09-25](https://www.checklistbank.org/dataset/303642/download)
    * Format: coldp
    * Choose root taxon: Fungi
    * Exclude ranks below: species
    * Include synonyms: Yes
    * Download and unzip this file (for example, your command might look like this, but your URL and filename will be different):
        ```
        wget https://download.checklistbank.org/job/01/01dc2cb7-3400-4f50-b36c-520236786a38.zip
        unzip 01dc2cb7-3400-4f50-b36c-520236786a38.zip
        ```
    * You should get two files (archived in this github repo as well). Subsequent steps:
        * [NameUsage.tsv](NameUsage.tsv)
        * [metadata.yaml](metadata.yaml)

2. Check how many different species name statuses there are in this file:

    ```
    tsvtk filter2 -f '$col:rank=="species"' NameUsage.tsv \
    | tsvtk cut -f "col:status" \
    | tsvtk freq
    ```

    In our case, we got the following:

    | col:status              | frequency|
    |-------------------------|--------|
    | provisionally accepted  |      28|
    | accepted                |  155841|
    | ambiguous synonym       |    3979|
    | synonym                 |  119404|

    i.e. there were 155,869 accepted and provisionally accepted species names, and 123,383 synonyms and ambiguous synonyms

3. Look for ECM genera names in both NameUsage.tsv and NameUsage_synonym_resolved.tsv:

    ```
    perl -plne 's/^/^/; s/$/\\s/' ECM_genera2.csv \
    | tsvtk grep -r -f 5 -P - \
        <(tsvtk filter2 -f '$col:rank=="species" && $col:status !~ "synonym"' NameUsage.tsv) -U \
    > ECM_accepted_match.tsv
    ```

    The `perl -plne 's/^/^/; s/$/\\s/'` bit ensures that the ECM genera names are turned into a regular expression which will only match a Genus name and not part of a genus name or part of a specific name.
    
    This list is then passed to the `tsvtk grep` command where `-r` specifies regular expressions on field `-f 5` and `-P -` says take the patterns from the stdin.