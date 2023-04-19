# GISAID ingestion for Flu PathMut

**NOTE: While this code is open-source, this workflow is intended for internal use only. It utilizes a GISAID data source that is not intended for use by the general public. Please see the `workflow_flu_genbank_ingest` to receive data from GenBank, or use either workflows as a template to ingest custom/in-house data.**

## Data Scraping

1. Navigate to GISAID EpiFlu portal
2. Select relevant serotype
   - Scraping has to be done per-serotype, so we can only select one at a time. For example, to get H1N1, select Type = A, H = 1, N = 1
3. Subset sequences by date
   - Can either do collection date or submission date, but submission date should be simpler and allow for less retroactive scraping (i.e., more likely for old isolates to be submitted, than for an old submission to be in limbo for a long time)
   - The maximum number of sequences that can be downloaded at one time is 10,000, so make sure that your selection yields <10,000 sequences. If you have more, then reduce your date range so that you get less than 10,000 sequences
   - Make a note of this date range—and put the date range into the name of your downloaded files so you can keep track of which ranges you have scraped so far.
   - Dates are inclusive, i.e., [start, end]
4. Once done subsetting, click "Search" at the bottom right
5. Select all sequences by clicking the checkbox at the top left of the table (this may take a few seconds)
6. Click "Download" at the bottom right - this should open up a dialog window
7. Select "Isolates as XLS (virus metadata only)" and click Download at the bottom right of the dialog.
   - This could take a while. Do not close the window or dialog box! (Still have to download sequences as well)
8. After the metadata file is downloaded, the dialog window should close. Rename the downloaded file to include the date range of this current selection
9. With the same sequences selected, click the "Download" button at the bottom right of the table again, and this time in the dialog window select "Sequences (DNA) as FASTA"
10. Select "All" for the "DNA" section (select all segments). The remaining default settings are ok.
11. Click "Download" at the bottom right of the dialog window.
    - Again, this could take a while
12. Rename the downloaded fasta file with the selected date range
13. Repeat the above steps, but for the next date range that will result in <10,000 sequences.

To speed up the process, you can open multiple tabs and run multiple downloads at once. Just make sure to not mix up your downloaded files!

14. Once all files are downloaded, move all metadata (.xls) files into the flu data `metadata/` folder and all sequence (.fasta) files into the flu data `sequences/` folder.

## Running

```
snakemake --cores 4
```

## Known Issues

If running into an error during the `copy_changed_files` step, try raising your open files limit with:

```
ulimit -n 100000
```
