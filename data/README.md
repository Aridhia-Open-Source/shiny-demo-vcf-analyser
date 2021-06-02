This folder comprises the data needed to run the VCF mini-apps.

It includes the VCF data however the clinvar data is a bit large at 100MB so this will need downloading and saving as e.g. clinvar_variant_summary.csv. It can
be downloaded (then unpacked) from the following FTP site at the time of writing:

ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

There are two VCF files included, however only one of these needs to be loaded (it depends how much data you want to use for your demo).

The VCF files were created from vcf_parser.py script in the top level 'genomics_visualisations' directory. This script splits VCF files into a header and a data CSV to
allow the data to be loaded into the platform and used in a relational manner. So it's quite a useful script but you may find it's still buggy while we
iron it out.