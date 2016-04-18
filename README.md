# Description

malign allows you to parse multiz alignment files (maf) and query the sequence in one species as aligned to a target species.  This project was born from a need to extract the sequence of mouse in the corresponding location in the human reference sequence in order to remove potential false positive sequence variants in xenograft models.  

Many people told me initially to use liftover [https://genome.ucsc.edu/cgi-bin/hgLiftOver] which does offer some between species conversions for example human to mouse.  However, in working with this approach I found a number of instances of no results being returned while I could tell by muliz alignment that the region was actually shared.  You might consider using liftover in place of this or in addition to it.

I could not find another utility to handle the maf files and do the conversions.  If a better tool exists please post a comment in issues and I'll take a look at it. 

# Requirements

You will have to have python2 installed, this is fairly standard on Linux platforms.  You will need to download the .maf.gz files from UCSC [http://hgdownload.cse.ucsc.edu/goldenPath/hg19/multiz100way/maf/];

```bash
rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/multiz100way/ ./
```

This doesn't appear to be versioned with this path, so you might wright down when you did the rsync or make sure that the file dates make sense.

If you go with the Jython option you should install that and Java as well.

# Usage

malign is pure Python and might work quite nicely in Jython if you plan to really hammer it with requests.  You start alignment_server.py with;

```bash
python alignment_server.py /path/to/maf/directory
```

The default port is 6480 and the default source and query are hg19 and mm10.  These can all be configured with command arguments.

Once you have the server started you can start to query it:

```bash
> curl -d '{"chrom": "M", "start": 5563, "end": 5564}' http://localhost:6480/sequence
{"action": "accepted", "queryseq": "CA", "sourceseq": "G-"}
```

The first time that a chromosome is queried it will load, which can take a very long time, subsequent queries are quite fast, so it wouldn't be such a crazy idea to just start one of these and have it running for some time if you had a server with loads of memory available if you planned to use this frequently.
