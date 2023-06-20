Download everything by clicking 'Clone or Download' > 'Download ZIP'.

See file 'Introduction_README.pdf' for instructions on using the application.  See below for notes on updates and issues.

6/20/2023-- NOTICE OF (MINOR) KNOWN ISSUES
There are four known issues that users of BPC should avoid when using it to analyze detrital geochronology datasets.  Note that these issues fall into one of two categories: Issues 1 and 2 are issues that will prevent the code from running due to formatting of data tables.  Issues 3 and 4 affect plotting BPC results using the built-in plotting tools.  The code also generates a spreadsheet table with BPC results, so BPC values are accessible even for users who experience issues 3 and 4.  It is important to note that none of these issues appear to influence the accuracy or reliability of BPC results.  I hope to make an update in the next few months to address these issues.  If you are interested in using the BPC method, I am happy to provide any needed support.  Please contact me at alex.tye[at]utahtech.edu

Issues:
1. Mac users can trigger an error depending on what type of comma separated value file format they use as inputs to the code.  When saving data tables for use with BPC, if using Excel, it is suggested that Mac users save the files as the "Comma Separated Values (Macintosh)" format rather than the "CSV UTF-8" format.  I have not encountered similar problems on Windows machines.
2. Uncertainty values of 0 will trigger an error.  Try changing any such values to small positive values.
3. Some users have had difficulty getting the function to change the order in which samples are plotted to work.  I have not encountered this problem but will look into it.  A potential workaround in the time being is to append a letter to the beginning of each sample name such that they will be automatically read in the desired order (e.g., "a14HSPK01, b15CRM03", etc...).
4. Some users have had difficulty getting the color scale to function properly for the output figure.  I will look into this.  In the meantime, users with programming experience can generate a table with color-coded values based on the table output values for BPC.
