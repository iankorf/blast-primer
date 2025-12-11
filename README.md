# BLAST Primer

## Contents

- Assumptions
- Setup
- Global and Local Alignment
- Random Expectation
- Local Alignment Statistics
- BLAST

## Assumptions

Sequence alignment exists at the intersection of genetics, molecular biology,
statistics, computer science, and information technology. If you're unfamiliar
with these topics, some parts of this will be confusing. This primer is meant
for my students and interns, so there isn't much hand-holding for the general
public. You will need to be comfortable with a Unix CLI to proceed.

This primer assumes you have an x86/AMD chipset running Linux or Windows or a
Mac. If you're using an unusual OS or hardware, some parts of the exercises may
fail.

## Setup

### Mac

Your computer is already running a Unix called Darwin. Open your `Terminal`
application to get to the CLI. If you haven't installed the command line tools
before, do so now.

```
xcode-select --install
```

### Windows

Either install the Windows Subsystem for Linux (WSL) via the Microsoft store
(it's free) or Install VirtualBox (also free). There are lots of Linux
distributions available either way. Popular choices include Ubuntu and Debian.
There are lots of guides and videos online.

### Install Conda

If you don't already have conda installed, install conda via the miniforge
installer. You can use `curl` or `wget`.

```
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
```

Don't do both of these.

```
wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
```

Run the shell script to install conda. Enter `yes` for the license agreement
and `yes` to install the base environment on startup.

```
bash Miniforge3-$(uname)-$(uname -m).sh
```

Open a new terminal. You should see `(base)` at the start of the CLI. If you do
not, seek help.


### Install EMBOSS and BLAST

The conda environment is specified in `blast-primer.yml`. Examine the file and
you will see that it specifies the environment name `blast-primer` and 2
packages: `emboss` and `blast-legacy`.

```
cat blast-primer.yml
```

Install the emboss and blast-legacy packages with the following command:

```
conda env create -f blast-primer.yml
```

Activate the `blast-primer` environment:

```
conda activate blast-primer
```

You now have access to several new programs in the `emboss` and `blast-legacy`
packages.


## Global and Local Alignment

There are two major classes of pairwise alignment: global and local. A global
alignment matches two sequences along their entire lengths while a local
alignment is patch of similarity. Let's observe each in action using the
`needle` and `water` programs from EMBOSS.

### Global Alignment Examples

Take a look at the simple sequences in the following FASTA files `s1.fa`,
`s2.fa`, and `s3.fa`.

```
cat s*.fa
>s1
WATER
>s2
WETTER
>s3
HALFWAYTHERE
```

Now let's align them with the `needle` program to observe how global alignment
works. After entering the command below, the program will ask you for gap
opening penalty and gap extension penalty. Accept its default values by hitting
the enter key.

```
needle -asequence s1.fa -bsequence s2.fa -out out
less out
```

Scroll down a bit to see the alignment. Every identical letter has a pipe
symbol between the two sequences, but the A and T have a `.` because they
mismatch. There is also a gap in sequence s1 to make up for the fact that the
two sequences aren't the same length.

```
s1                 1 W-ATER      5
                     | .|||
s2                 1 WETTER      6

```

Let's try another.

```
needle -asequence s1.fa -bsequence s3.fa -out out
less out
```

These two sequences are very different in length. Most of the alignment is
gaps. This is one reason global alignment isn't used very often.

```
s1                 1 ----WA--TER-      5
                         ||  .||
s3                 1 HALFWAYTHERE     12
```

### Local Alignment Examples

Let's try the same sequence alignments as before, but this time using the
`water` program which performs local alignment.

```
water -asequence s1.fa -bsequence s2.fa -out out
less out
```

In cases where the sequences are similar along their entire lengths, global and
local alignment can give the exact same alignment.

```
s1                 1 W-ATER      5
                     | .|||
s2                 1 WETTER      6
```

Now let's try the aligning sequences of very different lengths.

```
water -asequence s1.fa -bsequence s3.fa -out out
less out
```

As you can see, only letter aligned. Believe it or not, this is the best local
similarity between the two sequences, which is called the "maximum scoring
pair" or MSP. The `water` program is an implementation of the Smith-Waterman
algorithm, which finds a single MSP (even if there is more than one alignment
with the same score, the algorithm only returns one).

```
s1                 1 WA      2
                     ||
s3                 5 WA      6
```

You might be asking yourself why the local alignment isn't longer. Let's take a
closer look at the global alignment.

```
s1                 1 ----WA--TER-      5
                         ||  .||
s3                 1 HALFWAYTHERE     12
```

The gap characters result in a large score penalty. The only way to absorb the
penalty of a gap is to have enough matches on either side. The score of
aligning `WA` to `WA` is much greater than `TER` to `HER`. The `WA` alignment
gets 11 for the tryptophan matches and 4 for the Alanine matches for a total of
15 (which is the alignment score for `water`). The `TER` to `HER` match
generates -2, +5, +5. So in fact, it's better to consider `ER` as the matching
part. 15 is greater than 10 so that's why the best local alignment is just the
`WA`. If you want to see the scoring matrix used, you will find it here:

```
less ~/miniforge3/envs/blast-primer/data/BLOSUM62
```

You might be wondering why pairing `W` is worth more than pairing other amino
acids. It's because scoring matrix values are log-odds transformations of
observed divided by expected. `W` is a rare amino acid, so it's not expected to
occur very often. Therefore, when you see it, it's worth more.

## Random Expectation

Before doing an experiment, it's often useful to have a null hypothesis for
comparison. This is usually some form of random expectation. Let's look at what
random expectation looks like for sequence alignment, because it's probably not
what you expect. For this exercise, and the ones that follow, we will be using
nucleotide sequences with simple scoring schemes where all matches receive the
same score (unlike BLOSUM62 above).

Both global and local alignment are based on optimizing an alignment _score_.
You generally gain score when letters are the same and lose score when they are
different. In a simple scoring scheme, there are 3 types of scores.

- Match: reward for matching letters
- Mismatch: the penalty for mismatched letters
- Gap: the penalty for including a gap

In the global and local alignment examples above, we didn't pay much attention
to the alignment parameters (match, mismatch, gap), but it turns out that they
are ABSOLUTELY CRITICAL. Let's explore this with random sequences.

Take a look at the FASTA file `rseq.fa`. This contains two sequences that were
generated randomly using the script `randomseq` included in the repo.

```
>random-1
CACTATAGACATATAGGCTCACACGGGTTTGGGCTGCCGATATCCGGTAA
TCCTGCAAACTAGCACCTATCGACATGACAAACACCAGAATCTTACGCTG
>random-2
CTCCGGTGCCCCAAAGCTAACATTAGCGATATGGGATTTAGTGCGCCCAC
TTCATACCCCTCACTGCCAAATCTAAAGCAACGGTTACGCTCGGCGGGCC
```

Let's align these with different alignment parameters. First, let's try a very
simple scoring system: +1 match, -1 mismatch, -1 gap.

```
water -asequence rseq.fa -bsequence rseq.fa -gapopen 1 -gapextend 1 -datafile m+1-1.mat -out out
```

The `out` file contains two alignments, `random-1` aligned to `random-1` and
`random-1` aligned to `random-2`. It should be no surprise that when `random-1`
is aligned to itself, that it matches perfectly along its entire length with a
score of 100. But what about the other alignment 1 vs 2?

```
random-1           6 TAGACATATAGGCTCACACGGG-TTT-GGGCTGCCGA--T-AT-CCGGTA     49
                     || |||| || || .|.|.||| ||| |.|| |||.|  | || ||..|.
random-2          18 TA-ACAT-TA-GC-GATATGGGATTTAGTGC-GCCCACTTCATACCCCTC     62

random-1          50 ATCCTG-CAAA-CT--AGCA     65
                     |  ||| |||| ||  ||||
random-2          63 A--CTGCCAAATCTAAAGCA     80
```

This alignment is reported to have 64.3% identity with 45/70 matches and 17/70
gaps. Note that there is more than one way to calculate percent identity.

1. number of matches divided by length of alignment (count gaps)
2. number of matches divided by matches plus mismatches (ignore gaps)

Is the alignment above (1) 64.3% identical or (2) 84.9% identical? The problem
with (1) is that denominator can be larger than either sequence length. The
problem with (2) is that one shouldn't ignore gaps if they are present. Both
methods work similarly when the are few gaps, which is what happens most of the
time (provided you don't use idiotic alignment parameters). But it's worth
noting that something as simple as percent identity is up to interpretation.

### Varying Gap Cost

Keeping the match (+1) and mismatch (-1) scores the same, let's change the gap
costs and see what happens. Here is -1 gap (as shown previously). Note that I
use negative numbers to describe gaps even though the `water` program uses
positive values.

```
random-1           6 TAGACATATAGGCTCACACGGG-TTT-GGGCTGCCGA--T-AT-CCGGTA     49
                     || |||| || || .|.|.||| ||| |.|| |||.|  | || ||..|.
random-2          18 TA-ACAT-TA-GC-GATATGGGATTTAGTGC-GCCCACTTCATACCCCTC     62

random-1          50 ATCCTG-CAAA-CT--AGCA     65
                     |  ||| |||| ||  ||||
random-2          63 A--CTGCCAAATCTAAAGCA     80
```

Below is -2 gap. Note that the alignment is shorter, and there are no double
gaps anymore.

```
random-1          56 CAAA-CTAGCACCTATCGACATG     77
                     |||| |||.|| .||.|||.|||
random-2          12 CAAAGCTAACA-TTAGCGATATG     33
```

Below is -3 gap. This alignment is a little different from before, sacrificing
one of the gaps.

```
random-1          57 AAACTAGCACCTATCGACATG     77
                     ||.|||.|| .||.|||.|||
random-2          14 AAGCTAACA-TTAGCGATATG     33
```

Below is -4 gap or a penalty of greater magnitude. While the alignment below
doesn't contain any gaps, it still _could_. It would just need a bunch of
matches on either side to absorb the penalty of the gap.

```
random-1          93 TTACGCT     99
                     |||||||
random-2          85 TTACGCT     91
```

Which one of these alignments is the _correct_ one? They are all maximally
scoring under different scoring systems. When you choose a scoring system,
which you must do, you are affecting the alignment outcome.

### Affine Gaps

The previous examples used the same gap open and extension cost. That is a gap
that is 10 nt long costs 10 times that of a single gap. Graphing constant gap
costs result in a linear relationship going through the origin. In other words:
y = mx. Biological sequence gaps tend to show up in clusters, so the cost of 4
gaps might not be much more than 3. In sequence alignment, this is modeled with
_affine_ gap costs. In other words: y = mx + b. Let's see what happens when the
gap costs don't go through the origin. Here is what happens with gapopen 2 and
gapextend 1. As you can see, it's not the same as any of the previous
alignments.

```
random-1          43 TCCGGTAATCCTGCAAA-CTAGCACCTATCGACATG     77
                     ||||||  .|| .|||| |||.|| .||.|||.|||
random-2           2 TCCGGT--GCC-CCAAAGCTAACA-TTAGCGATATG     33
```

So is it better to have constant or affine gaps? Again, you get to choose. But
what should you choose?

### Ungapped Alignment

The examples above kept the match (+1) and mismatch (-1) scores the same and
changed the gap costs. Let's see what happens if we change the match and
mismatch scores while leaving the gap costs at -10. With such high gap
penalties, the there won't be any gaps. This will give us some insight into the
ungapped similarities of random sequences. This isn't just an exercise. It
turns out, that the entire foundation of local alignment statistics is based on
random ungapped alignments.

The `evd` program in this repo creates random sequences and aligns them with
`water`. Here's an example command line that creates 2 random sequences of
length 2000, and then aligns them with a +1/-1 scoring scheme.

```
./evd --length 2000 --matrix m+1-1.mat --count 2 --verbose
```

Each time you run the program it creates different random sequences, so don't
expect to see the alignment below. The program also reports the score, length,
and percent identity of the alignment.

```
r0              1974 ATAAAGGACAGATATGGCCGGGCA   1997
                     |||.|.|||||.|||.|||||.||
r1               553 ATATATGACAGTTATTGCCGGACA    576
```

Try running this a few times and you'll see a bunch of different alignments.
Recall that the whole point of this is to see what random similarity looks
like. In order to understand this random background better, let's perform this
experiment 1000 times. This should take about a minute.

```
./evd --length 2000 --matrix m+1-1.mat --count 1000 > out
```

The file `out` contains all of the scores, lengths, and percent identities of
the alignments. While you might expect the scores to follow a normal
distribution, they do not. Alignment is a maximizing operation, so the scores
are not bell curve, but a smushed bell curve called an extreme value
distribution (EVD). Let's make a quick histogram to observe.

```
cut -f1 out | sort -n | uniq -c
```

As you can see, the mode is at 12 and there is a long tail. Even though most
alignments have scores from 12.0 to 14.0, it's possible to get really high
scores like 20. However it's very difficult for the maximum score to be a low
score, like 9.

```
     85 11.0
    363 12.0
    315 13.0
    151 14.0
     50 15.0
     20 16.0
     12 17.0
      1 19.0
      2 20.0
```

Let's examine the average score, average length, and average percent identity.

```
awk '{t1+=$1; t2+=$2; t3+=$3; n++} END {print t1/n, t2/n, t3/n }' out
```

That `awk` script does some quick math. Don't worry if `awk` is unfamiliar to
you. You could do this in any programming language quickly enough. Here's the
answer I got (which will be a little different from you).

```
12.8458 22.8659 80.2933
```

The average alignment has a score of 12.8, a length of 22.9, and about 80%
identity. Now comes the really interesting part. Let's change the mismatch
parameter and examine the properties of the alignments. The table below
summarizes a couple runs of `evd` with different match and mismatch values.

| Match | Mismatch | Score | Length | Percent |
| :---: | :------: | :---: | :----: | :-----: |
|  +1   |    -1    | 11.57 | 20.08  |  81.31  |
|  +1   |    -2    | 9.93  | 11.14  |  97.07  |
|  +1   |    -3    | 9.72  | 10.00  |  99.48  |

**The choice of match and mismatch determines what you find.**

If you're looking for alignments that are 80% identical, +1/-1 is a good
choice. If you're looking for 99% identity, +1/-3 is a good choice. If you use
the default values, you might not get the behavior you want. Every alignment
program has its own default values. BLASTN has used +5/-4, +1/-3, and +1/-2
over the course of its development, so depending on what version you use, you
may get very different default behaviors.

### Gapped Alignment

Let's do the same `evd` experiments, but this time we'll add various gap
opening and extension values. Here, we'll keep the +1/-2 scoring system
constant.

| Match | Mismatch | GapO | GapE | Score | Length | Percent |
| :---: | :------: | :--: | :--: | :---: | :----: | :-----: |
|  +1   |    -2    |  -3  |  -3  | 10.06 | 12.13  |  95.87  |
|  +1   |    -2    |  -2  |  -2  | 10.65 | 16.20  |  90.52  |
|  +1   |    -2    |  -2  |  -1  | 10.90 | 19.12  |  87.53  |
|  +1   |    -1    |  -2  |  -2  | 13.97 | 49.11  |  70.97  |

Notice that the smaller the gap penalties, the longer and lower percent
identity of the alignment. It's easier to make alignments if you allow gaps,
and the lower the gap costs, the easier it gets.

### Synthesis

Every scoring system has a _stringency_. Some scoring systems are stringent
while others are permissive. One way to describe the stringency is by the
average percent identity of random alignments.

| Match | Mismatch | GapO | GapE | Pct  | Stringency |
| :---: | :------: | :--: | :--: | :--: | :--------- |
|  +1   |    -3    |      |      | 99.5 | very high  |
|  +1   |    -2    |      |      | 97.1 | high       |
|  +1   |    -2    |  -3  |  -3  | 95.9 | high       |
|  +1   |    -2    |  -2  |  -2  | 90.5 | moderate   |
|  +1   |    -2    |  -2  |  -1  | 87.5 | moderate   |
|  +1   |    -1    |      |      | 81.3 | low        |
|  +1   |    -1    |  -2  |  -2  | 71.0 | very low   |

The original version of BLAST defaulted to +5/-4 with no gapping allowed.
Later, this was +1/-3 with gaps -5/-2. Later, it was +1/-2 with gaps -2.5/-2.5.
What is it now? You'll have to check the command line options. EMBOSS `water`
defaults to +5/-4 with gaps -10/-0.5. Clearly, despite how critical they are,
there is no consensus on the proper default alignment parameters. This is sort
of like hybridization kinetics. There is no one setting for a PCR machine that
works for all oligos. Each person must set up the machine for their primers.
And so it should be with BLAST.

Note that there are 2 ways to describe gap costs. NCBI-BLAST versions 2.0+
follow a linear formula as y = mx + b, where y is the total gap cost, b is the
gap opening cost, and m is the cost of extending gaps. AB-BLAST (and its
historical descendants WU-BLAST and NCBI-BLAST 1.4) and EMBOSS `water` and
`needle` describe the cost of the first gap followed by the cost of additional
gaps. So, NCBI 3/1 is the same as `water` 4/1.

### Summary

The choice of alignment parameters (match, mismatch, gap) has a major effect on
the random background of pairwise alignment. Performing sequence alignment
experiments without considering alignment parameters is sort of like baking and
not caring about the temperature of the oven.

## Local Alignment Statistics

Any two sequences compared by Smith-Waterman will have a maximum scoring pair.
It doesn't matter if the sequences are related or not: there will always be an
MSP. By analogy, you could compare two completely different books, a fantasy
novel written in English and a home repair guide in French, and they would
still have an MSP. So the next logical question is this: how can you tell if
the MSP is any good?

One of the important innovations of BLAST is that it implemented
Karlin-Altschul statistics. The K-A equation tells us how often a score can
happen by chance. For example, it can tell you that a score of 10 is supposed
to occur very often and a score of 100 is not. Here's the equation.

> E = kMNe^-ls

- E: the expected number of alignments by chance (e-value)
- k: a minor constant we won't discuss
- M: the size of one sequence
- N: the size of the other sequence
- e: 2.71828...
- l: lambda is the scaling factor of the matrix
- s: the score of the alignment

To get more intuition on how the K-A equation works, let's simplify it by
dropping k and lambda, and rounding e to 2.

> E = MN2^-s

The number of alignments expected by chance depends on the search space (MN)
and the score of an alignment (s). If we double the search space (e.g. 2M), we
double the number of alignments expected by chance. If we decrease the score of
an alignment by 1, we also double the number of alignments by chance. We can
gain several really important intuitions from this:

- Low scoring alignments can occur frequently
- High scoring alignments do not occur by chance
- Huge search spaces may contain random, high scoring alignments

The K-A equation makes a few assumptions:

1. A positive score must be possible
2. The expected score must be negative
3. The letters in the sequences are independent and identically distributed
4. The sequences are infinitely long
5. The alignments do not contain gaps

Many of these assumptions are violated by real biological sequences. Of course
no sequences are infinitely long. One of the most troublesome assumptions is
that alignments have no gaps. The reason BLAST was originally an ungapped
alignment algorithm was because it used the K-A equation to estimate the
significance of alignment scores.

Lambda is a critical parameter that normalizes different scoring schemes. You
can imagine that a +1/-1 scoring scheme should behave very similarly to +2/-2.
In fact, they behave identically. The difference is that +1/-1 has a lambda
twice as large as +2/-2. As a result, it doesn't matter if you use a +1/-1 or
+2/-2, the e-values will be exactly the same.

It turns out that the K-A equation can be applied to gapped alignments
by borrowing the lambda from an equally stringent scoring scheme.
Basically, you do a bunch of random gapped alignments and look for
similar percent identities to ungapped alignments. However, BLAST
doesn't support all possible gapped scoring schemes. You have to use a
scoring scheme it already knows about.

## BLAST

### bl2seq

Bl2seq performs a comparison between two sequences (either proteins or
nucleotides) using either the blastn or blastp algorithm. The command
compares a sequence against either a local database or a second sequence.

Let's compare the two actin proteins from Gallus gallus (chicken) and
Drosophila melanogaster (fruit fly). They are in the files `dm.fa` and `gg.fa`.

```
bl2seq -i gg.fa -j dm.fa -p blastp
```

Once this command is run, you should see an output detailing the alignment of the two sequences.

```
Query: 2   DDDIAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSK 61
           D+++AALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSK
Sbjct: 3   DEEVAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSK 62
```

```
Query: 122 IMFETFNTPAMYVAIQAVLSLYASGRTTGIVMDSGDGVTHTVPIYEGYA----LPHAILR 177
           IMFETFNTPAMYVAIQAVLSLYASGRTTGIV+DSGDGV+HTVPIYEGYA    LPHAILR
Sbjct: 120 IMFETFNTPAMYVAIQAVLSLYASGRTTGIVLDSGDGVSHTVPIYEGYAAAAALPHAILR 179
```

In these two alignments, pluses, minuses, and blank spaces are used to
represent different aspects of the alignment between the query and the
subject.

The plus (+) symbol indicates positions where the amino acids in the
aligned sequences are similar. The minus (-) symbol indicates positions
where there is an insertion or deletion in one of the sequences. A blank
space indicates positions in the alignnment where the amino acids have no
match or similarities.

Bl2seq also outputs Lambda, K, and H (variables used to calculate the
statistics of alignment). These variables are used to determine the
significance of the alignment which helps us understand the validity of its
MSP.

In the alignment above the Lambda, K, and H variables should be the
following:

```
Lambda     K      H
   0.319    0.135    0.401

Gapped
Lambda     K      H
   0.267   0.0410    0.140
```

### blast-legacy

BLAST is usually used to search/align one (or more) sequences to a databases of
sequences. For example, we could search all of the E.coli proteins vs. all of
the Y.pestis proteins (Yersinia pestis is the bacterium responsible for bubonic
plague). To do that, one of the fasta files must be turned into the database.
The command that turns fasta files into blast-able databases is called
`formatdb`.

```
formatdb -i E.coli.faa
```

After running `formatdb`, you should see some new files with extensions
`.phr`, `.psq`, and `.pin`. These files comprise the blast database.

Our stated goal was to align all of the Y.pestis proteins to E.coli.
However, we don't have any idea how long that will take. It could take
seconds or hours. So the first thing we should do is make a subset of
the Y.pestis proteome and search that against E.coli. A simple `wc` will tell us
how many lines are in the file.

```
wc Y.pestis.faa
```

There are about 24k lines. Let's make a smaller version of Y.pestis containing
about 1% of that total.

```
head -2400 Y.pestis.faa > mini.faa
```

A simple `grep` verifies that there are around 1% of the proteins in the
`mini.faa` compared to the full proteome.

```
grep -c ">" Y.pestis.faa mini.faa
```

Now it's finally time to run BLAST. Since both the query and database
are proteins, the program type is `blastp`. We'll save the output in a
throwaway file called `foo`. We'll prepend the command with `time` to
monitor resources. Here's the command line:

```
time blastall -p blastp -d E.coli.faa -i mini.faa > foo
```

Given how long it took for `water` to align 1000 sequences to each other, you
might think it will take `blastp` a long time to align 494 Y.pestis proteins to
4298 E.coli proteins. It does't. This should take under 10 seconds (depending
on your hardware). You will see a few warnings about Selenocysteine that you
can ignore. Examine the output file with `less foo`. By default, BLAST reports
alignment with high E-values, which represent random similarities. In the
future, we will set `-e 1e-5` to remove the really poor alignments.

By default `blastall` uses only 1 cpu. You can speed up the search by
giving it more CPUs with the `-a` parameter. Here are the results I got
with additional CPUs.

| CPUs | Time | CLI
|:----:|:----:|:-------------------------------------------------------
|   1  | 2.17 | `blastp -d E.coli.faa -i mini.faa -e 1e-5 -a 1 > foo`
|   2  | 1.31 | `blastp -d E.coli.faa -i mini.faa -e 1e-5 -a 2 > foo`
|   3  | 1.08 | `blastp -d E.coli.faa -i mini.faa -e 1e-5 -a 3 > foo`
|   4  | 0.92 | `blastp -d E.coli.faa -i mini.faa -e 1e-5 -a 4 > foo`

As you can see, there are diminishing returns with more CPUs. Some parts
of BLAST cannot be parallelized across multiple CPUs. For example, all
of the CPUs must eventually write their output, and this is not
parallelized.

Now it's time to search the whole Y.pestis proteome against the whole
E.coli proteome. We can now estimate how long the job will take: about
100x longer with the full Y.pestis proteome. In reality, it may take
longer or shorter than the estimate, but at least you have some idea.

```
blastall -p blastp -d E.coli.faa -i Y.pestis.faa -e 1e-5 -a 4 > yve.blastp
```

### Parsing a BLAST report

The `yve.blastp` output file is 28M. That's a lot of text to look through for a
human. If you don't need to examine the alignments, you can get the numerical
data by changing the output format to tabular using `-m 9`. This reduces the
file to just 3M. Let's try that. Also we'll add another parameter `-f 13` to
make the search go faster at the cost of sensitivity.

```
blastall -p blastp -d E.coli.faa -i Y.pestis.faa -e 1e-5 -a 4 -m 9 -f 13 > yve.tsv
```

Page through the outputs with `less` to get a feeling for what a proteome to
proteome alignment looks like.

```
less yve.blastp
less yve.tsv
```

You might notice that some query sequences from Y.pestis have multiple matches
in E.coli and some have none. How could you display the number of matches?

```
grep -v "^#" yve.tsv | less | cut -f1 | uniq -c | less
```

Note that this doesn't show any zeroes. If a Y.pestis protein doesn't match
anything in E.coli, it's not in the blast report. What if the genes that make
Y.pestis "bad" are genes that are specifically not in E.coli?

Look up the `comm` command and then use that to determine how many proteins in
Y.pestis are not in E.coli.
