# VarGenius - Frequently Asked Questions


Here you can find issues already solved at any step with *VarGenius*. Please search here before to post any question on the forum.
   

*Question*: *VarGenius* gave the following error:
```
Exception in thread "main" java.lang.UnsupportedClassVersionError: org/broadinstitute/gatk/engine/CommandLineGATK : Unsupported major.minor version 52.0
```

*Answer*: This is an error that happens with GATK when using an uncompatible java version. Please search the GATK website the best java version to use for GATK.


-------------------------

*Question*: *VarGenius* gave the following error:
```
Error: the required database file hg19_cadd13.txt does not exist
```

*Answer*: This is an error occurring when Annovar cannot download a database. You can try to re-run *VarGenius* using **annov_db_dl = YES** or try to install manually the database using **annotate_variation.pl -downdb  -webfrom annovar -buildver hg19 cadd13  /..pathto../annovar/humandb**

-------------------------


------------------------------------------------

If you get some error during the installation or the running of *VarGenius* please see the FAQ page of the user guide or ask on the google group: https://groups.google.com/forum/#!forum/VarGenius
