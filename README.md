# VCF Analyser

Variant Calling Format (VCF) files are used for storing gene sequence variations. Much of the gneomic data can be redundant, as the genetic information is shared accross multiple genomes, VCF files only show the variations along with a reference genome. VCF are composed by a header and a body, the header provides metadata describing the rest of the file. The body is tab separated into 8 mandatory columns plus a number of other optional columns. The mandatory columns are: 

- CHROM - Chromosome number
- POS - Position of the variation
- ID - Identifier of the variation, if unknown ".
- REF - Base in the reference sequence
- ALT - Alternative bases in the position
- QUAL - Quality score
- FILTER - What filters the variation has passed
- INFO - description of the variation

This mini-app offers a basic variant exploratory data analysis.

- Exploratory analysis of variants
- Visualize patterns (mutation spread, quality scores)
- Link variants to a data source
- Basic QC analysis

### Variant Browser

The table in the middle of the page shows all the body of the VCF file. This table can be filtered by sample, chromosome and position using the options in the left-side column (Inputs). If you select one variant by clicking on the row, the right-side box be filled with the annotation information for that variant. The information displayed consist of:

- Coordinates (CHR:POS)
- Link to open the variant in Ensemble
- Base change between Reference (REF) and Alternative (ALT)
- Quality
- Number of clinical variants

A small table will also appear with more information about the selected variant.

If you click on "View VCF Header", the header of the VCF file will appear at the bottom of the screen.

### SNP Densities

This tab displays a graph showing the SNP distributions. Using the left-side box (Inputs), the user can filter the data being displayed in the graph by choosing a different chromosome and a specific genomic range.

### SNP Type Distribution

This tab displays a heatmap showing the number of SNPs by base change. Using the left-side box (Inputs), the user can filter the data being displayed in the graph by choosing a different chromosome and a specific genomic range.


## Data

The data for this demo app was obtained using the FTP site in ClinVar (https://www.ncbi.nlm.nih.gov/clinvar/). At the time of writing: ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz.

## Checkout and Run

You can clone this repository by using the command:

```
git clone https://github.com/aridhia/demo-vcf/analyser
```

Open the .Rproj file in RStudio, soruce the script called `dependencies.R` to install all the packages required by the app and run `runApp()` to start the app.

### Deploying to the workspace

1. Download this GitHub repo as a .zip file.
2. Create a new blank Shiny app in your workspace called "vcf-analyser".
3. Navigate to the `vcf-analyser` folder under "files".
4. Delete the `app.R` file from the `vcf-analyser` folder. Make sure you keep the `.version` file!
5. Upload the .zip file to the `vcf-analyser` folder.
6. Extract the .zip file. Make sure "Folder name" is blank and "Remove compressed file after extracting" is ticked.
7. Navigate into the unzipped folder.
8. Select all content of the unzipped folder, and move it to the `vcf-analyser` folder (so, one level up).
9. Delete the now empty unzipped folder.
10. Start the R console and run the `dependencies.R` script to install all R packages that the app requires.
11. Run the app in your workspace.

For more information visit https://knowledgebase.aridhia.io/article/how-to-upload-your-mini-app/
