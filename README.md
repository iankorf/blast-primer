BLAST Primer
============

## Contents ##

+ Assumptions
+ Setup
+ Global and Local Alignment
+ Random Expectation
+ Local Alignment Statistics
+ BLAST

## Assumptions ##

Sequence alignment exists at the intersection of genetics, molecular biology,
statistics, computer science, and information technology. If you're unfamiliar
with these topics, some parts of this will be confusing. This primer is meant
for my students and interns, so there isn't much hand-holding for the general
public. You will need to be comfortable with a Unix CLI to proceed.

This primer assumes you have an x86/AMD chipset running Linux (either native or
virtualized in Windows) or an older Mac with an Intel chipset. If any of the
following situations apply to you, expect some pain if not outright failure.

+ Newer (post-2020-ish) Mac with Apple Silicon chips
+ Windows running Windows Substem for Linux (WSL)
+ Windows running Cygwin
+ Windows running GitBash
+ ChromeOS
+ Raspberry Pi or other SBC

## Setup ##

### Intel-Mac

If you have an older Mac (pre-2020-ish), you already have a Unix CLI. Just open
your `Terminal` application.

### Windows

+ Install VirtualBox
+ Create a VM with 2-4G RAM and a flexible size drive of 40G
+ Install a lightweight Linux distribution like Lubuntu, LinuxLite, or Mint

### Miniconda

Install Miniconda and then perform the following post-installtion tweaks to
ensure packages from BioConda are installed correctly.

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

### Install EMBOSS and BLAST

Create a conda environment for this primer and install emboss and 2 different
versions of blast. This will take a couple minutes to install, and will end up
consuming about 3G of space.

```
conda create --name blast-primer emboss blast-legacy blast
conda activate blast-primer
```

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
THEWEATHERISCALM
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
s1                 1 ---W-ATER-------      5
                        | ||..
s3                 1 THEWEATHERISCALM     16
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
algorith, which finds a single MSP (even if there is more than one alignment
with the same score, the algorithm only returns one).

```
s1                 1 W      1
                     |
s3                 4 W      4
```

You might be asking yourself why the local alignment isn't longer. It's because
of the default gap penalties we accepted. You might also be asking why aligning
a single `W` in each sequence is better than aligning the pair of `AT`s found
in both `WATER` and `WEATHER`. It's because the score of alinging tryptophan
(W) is greater than aligning both an alanine (A) and threonine (T). The reason
why W has a greater score is because scores are essentially log odds ratios of
observed / expected and tryptophan is a very rare amino acid.

## Random Expectation

Before doing an experiment, it's often useful to have a null hypothesis for
comparison. This is usually some form of random expectation. Let's look at what
random expectation looks like for sequence alignment, because it's probably not
what you expect. For this exercise, and the ones that follow, we will be using
nucleotide sequences with simple scoring schemes where all matches receive the
same score (unlike the tryptophan story above).

Both global and local alignment are based on optimizing an alignent _score_.
You generally gain score when letters are the same and lose score when they are
different. In a simple scoring scheme, there are 3 types of scores.

+ Match: reward for matching letters
+ Mismatch: the penalty for mismatched letters
+ Gap: the penalty for including a gap

In the global and local alignment examples above, we didn't pay much attention
to the alignment parameters (match, mismatch, gap), but it turns out that they
are ABSOLUTELY CRITICAL. Let's explore with random sequences.

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
score of 100. But what about the other?

```
random-1           6 TAGACATATAGGCTCACACGGG-TTT-GGGCTGCCGA--T-AT-CCGGTA     49
                     || |||| || || .|.|.||| ||| |.|| |||.|  | || ||..|.
random-2          18 TA-ACAT-TA-GC-GATATGGGATTTAGTGC-GCCCACTTCATACCCCTC     62

random-1          50 ATCCTG-CAAA-CT--AGCA     65
                     |  ||| |||| ||  ||||
random-2          63 A--CTGCCAAATCTAAAGCA     80
```

This alignment is reported to have 64.3% identity with 45/70 matches and 17/70
gaps. There is more than one way to calculate percent identity.

1. number of matches divided by length of alignment (count gaps)
2. number of matches divided by matches plus mismatches (ignore gaps)

Is the alignment above (1) 64.3% identical or (2) 84.9% identical? The problem
with (1) is that denominator can be larger than either sequence length. The
problem with (2) is that one shouldn't ignore gaps if they are present. Both
methods work similarly when the are few gaps, which is what happens most of the
time (provided you don't use idiotic alignment parameters).

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

### Complex Gaps

The previous examples used the same gap open and extension cost. That is a gap
that is 10 nt long costs 10 times that of a single gap. Sequence gaps tend to
come in clusters, so we often use a two gap costs, one to open, and one to
extend. Let's see what happens when we vary the two.

Here is  what happens with gapopen 2 and gapextend 1. As you can see, it's
not the same as any of the previous alignments.

```
random-1          43 TCCGGTAATCCTGCAAA-CTAGCACCTATCGACATG     77
                     ||||||  .|| .|||| |||.|| .||.|||.|||
random-2           2 TCCGGT--GCC-CCAAAGCTAACA-TTAGCGATATG     33
```

### Ungapped Alignment

The examples above kept the match (+1) and mismatch (-1) scores the same and
changed the gap costs. Let's see what happens if we change the match and
mismatch scores while leaving the gap costs at -10. With such high gap
penalties, the there won't be many gaps. This will give us some insight into
the ungapped similarities of random sequences. This isn't just an exercise. It
turns out, that the entire foundation of local alignment statistics is based on
random ungapped alignments.

The `evd` program in this repo creates random sequences and aligns them with
`water`. Here's an example command line that creates 2 random sequences of
length 2000, aligns them with a +1/-1 scoring scheme.

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
experiment 1000 times.

```
/evd --length 2000 --matrix m+1-1.mat --count 1000 > out
```

The file `out` contains all of the scores, lengths, and percent identities of
the alignments. While you might expect the scores to follow a normal
distribution, they do not. Alignment is a maximizing operation, so the scores
are not bell curve, but a smushed bell curve called an extreme value
distribution (EVD). Let's make a quick histogram to observe.

```
cut -f out | sort -n | uniq -c
```

As you can see, the mode is at 12 and there is a long tail. Even though most
alignments have scores from 12.0 to 14.0, it's possible to get really high
scores like 20. However it's very difficult for the maximum score to be a very
low score.

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
12.8458 22.8659 80.2933
```

The average alignment has a score of 12.8, a length of 22.9, and about 80%
identity. Now comes the really interesting part. Let's change the alignment
parameters and examine the properties of the alignments. The table below
summarizes a couple runs of `evd` with different match and mismatch values.

| Match | Mismatch | Score | Length | Percent |
|:-----:|:--------:|:-----:|:------:|:-------:|
|   +1  |    -1    | 11.57 |  20.08 |  81.31  |
|   +1  |    -2    |  9.93 |  11.14 |  97.07  |
|   +1  |    -3    |  9.72 |  10.00 |  99.48  |

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
|:-----:|:--------:|:----:|:----:|:-----:|:------:|:-------:|
|   +1  |    -2    |  -3  |  -3  | 10.06 |  12.13 |  95.87  |
|   +1  |    -2    |  -2  |  -2  | 10.65 |  16.20 |  90.52  |
|   +1  |    -2    |  -2  |  -1  | 10.90 |  19.12 |  87.53  |
|   +1  |    -1    |  -2  |  -2  | 13.97 |  49.11 |  70.97  |

Notice that the smaller the gap penalties, the longer and lower pecent identity
of the alignment.

### Synthesis

Every scoring system has a _stringency_. Some scoring systems are stringent
while others are permissive. One way to describe the stringency is by the
average percent identity of random alignments.

| Match | Mismatch | GapO | GapE |  Pct  | Stringency
|:-----:|:--------:|:----:|:----:|:-----:|:-----------
|   +1  |    -3    |      |      |  99.5 | very high
|   +1  |    -2    |      |      |  97.1 | high
|   +1  |    -2    |  -3  |  -3  |  95.9 | high
|   +1  |    -2    |  -2  |  -2  |  90.5 | moderate
|   +1  |    -2    |  -2  |  -1  |  87.5 | moderate
|   +1  |    -1    |      |      |  81.3 | low
|   +1  |    -1    |  -2  |  -2  |  71.0 | very low

The orignal version of BLAST defaulted to +5/-4 with no gapping allowed. Later,
this was +1/-3 with gaps -5/-2. The current version is +1/-2 with gaps
-2.5/-2.5. EMBOSS `water` defaults to +5/-4 with gaps -10/-0.5. Clearly,
despite how critical they are, there is no consensus on the proper default
alignment parameters.

Note that there are 2 ways to describe gap costs. NCBI-BLAST versions 2.0+
follow a linear formula as y = mx + b, where y is the total gap cost, b is the
gap opening cost, and m is the cost of extending gaps. AB-BLAST (and its
historical descendents WU-BLAST and NCBI-BLAST 1.4) and EMBOSS `water` and
`needle` describe the cost of the first gap followed by the cost of additional
gaps. So the current NCBI-BLASTN scoring scheme is 0/2.5 in their parlance, but
2.5/2.5 in AB-BLAST and EMBOSS.

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
still have an MSP. So the next logical question is this: is the MSP any good?

One of the important innovations of BLAST is that it implemented
Karlin-Altschul statistics. The K-A equation tells us how often a score can
happen by chance. For example, it can tell you that a score of 10 is supposed
to occur very often and a score of 100 is not. Here's the equation.

>E = kMNe^-ls

+ E: the expected number of alignments by chance (e-value)
+ k: a minor constant we won't discuss
+ M: the size of one sequence
+ N: the size of the other sequence
+ e: 2.71828...
+ l: lambda is the scaling factor of the matrix
+ s: the score of the alignment

To get more intuition on how the K-A equation works, let's simplify it by
dropping k and lambda, and rounding e to 2.

>E = MN2^-s

The number of alignments expected by chance depends on the search space (MN)
and the score of an alignment (s). If we double the search space (e.g. 2M), we
double the number of alignments expected by chance. If we decrease the score of
an alignment by 1, we also double the number of alignments by chance. We can
gain several really important intuitions from this:

+ Low scoring alignments can occur frequently
+ High scoring alignments do not occur by chance
+ Huge search spaces may contain random, high scoring alignments

The K-A equation makes a few assumptions:

1. A positive score must be possible
2. The expected score must be negative
3. The letters in the sequences are independent and identially distributed
4. The sequences are infinitely long
5. The alignments do not contain gaps

Many of these assumptions are violoated by real biological sequences. Of course
no sequences are infinitely long. One of the most troublesome assumptions is
that alignments have no gaps. The reason BLAST was originally an ungapped
alignment algorithm was because it used the K-A equation to estimate the
significance of alignment scores.

Lambda is a critical parameter that normalizes different scoring schemes. You
can imagine that a +1/-1 scoring scheme should behave very similarly to +2/-2.
In fact, they behave identically. The difference is that +1/-1 has a lambda
twice as large as +2/-2. As a result, it doesn't matter if you use a +1/-1 or
+2/-2, the e-values will be exactly the same.

It turns out that the K-A equation can be applied to gapped alignments by
borrowing the lambda from an equally stringent scoring scheme. Basically, you
do a bunch of random gapped alignments and look for similar percent identities
to ungapped alignments.

## BLAST

+ multiple HSPs
+ seeding
+ extension
+ masking
	+ complexity filters
	+ hard
	+ soft
	+ word
+ gapped vs ungapped
