#Author: Shivam Patel

from typing import Any
import pandas as pd
import streamlit as st
import altair as alt
from PIL import Image


# Title and logo
image = Image.open('image.png')
st.image(image, use_column_width = True)

# Mission and Purpose
st.write("""

***

# Mission & Purpose

Dna Hub is an web-based application designed to analyze DNA for either clinical research or commercial product design.
User no longer needs to go through the tedious process of traveling to various sites to utilize essential tools. DNA Hub
is designed for have a central location for essential tools whether it be transcription, translation, nucleotide distrubtion, amino acid
distribution, or data on basic structure and interaction. DNA Hub aims to make bioinformatic research efficient! All you need is a FASTA Sequence.

****
""")

#BOX TO INPUT SEQUENCE
st.header("Enter desired DNA sequence you wish to perform analysis.")
sequence_input = ">DNA Query\nGTCGATGCATGTCGATGCATGCTAGTACGTAGTGACTGACTAGTCGATGCATGTCGATGCATGCTAGTACGTAGTGACTGACTAGTCGATGCATGTCGATGCATGCTAGTACGTAGTGACTGACTAGTCGATGCATGTCGATGCATGCTAGTACGTAGTGACTGACTAGTCGATGCATGTCGATGCATGCTAGTACGTAGTGACTGACTAGTCGATGCATGTCGATGCATGCTAGTACGTAGTGACTGACTAGTCGATGCATGTCGATGCATGCTAGTACGTAGTGACTGACTA"
sequence = st.text_area("Only FASTA file format is accepted. See sample sequence below. ", sequence_input, height=400)
#creates a list of each line within the query
sequence = sequence.splitlines()
#skips the first sequence in the list
sequence = sequence[1:]
sequence = ''.join(sequence)

st.write("""
***
""")

# INPUT DISPLAY
st.header("INPUT")
sequence

st.write("""
***
""")

# OUTPUT  DISPLAY
st.header("OUTPUT")

# VALIDATION
def validate(seq1):
    for x in seq1:
        test = False
        if x != 'A' and x != 'G' and x != 'T' and x != 'C' and x != '':
            st.write('This is not a valid sequence, please check sequence again')
            test = False
            break
        else:
            test = True
    return test


if validate(sequence) == True:
    st.write(sequence)
else:
    pass



# Nucleotide count
def nucleotide_count(seq):
    d = dict([
        ('A', seq.count('A')),
        ('T', seq.count('T')),
        ('G', seq.count('G')),
        ('C', seq.count('C'))

    ])
    return d
X = nucleotide_count(sequence)

st.subheader("Nucleotide Count")

#DATAFRAME
df = pd.DataFrame.from_dict(X, orient='index')
df = df.rename({0: 'count'}, axis ='columns')
df.reset_index(inplace=True)
df = df.rename(columns = {'index':'nucleotide'})
st.write(df)

#NUCELOTIDE DISTRIBUTION BAR PLOT

st.write("""

""")
p = alt.Chart(df).mark_bar().encode(
    x = 'nucleotide',
    y ='count').configure_mark(opacity = 0.2, color = 'gold')
p = p.properties(width= alt.Step(100))

st.write(p)
st.write("""
***
""")

# HYDROGEN BONDING
st.subheader("Hydrogen Bonding")

GC = (X['G'] + X['C'])
AT = (X['A'] + X['T'])
totalnucleotide = GC + AT
bonds = (GC * 3) + (AT *2)
GCratio = (GC/totalnucleotide) * 100
st.write(' The GC content indicates ' + str(GCratio) + '%. ' +
         'The GC content refrences structural stability of DNA. ' +
         'There is ' + str(bonds) + ' hydrogen bonds in this DNA. ' +
         'This includes ' + str(GC) + ' triple bonds from guanine and cytosine and ' +
         str(AT) + ' double bonds from adenine and thymine.'


         )



def stability(ratio):
    if ratio < 40:
        st.write('GC content of this sequence is too low, and potentially unstable DNA for Primer to bind.')
    elif ratio >60:
        st.write('GC content of this sequence is too high, and can \n ' +
                 'affect the efficiency of PCR due to the tendency of these templates to fold ' +
                 'into complex secondary structures. \n It will directly affect the melting ' +
                 'temperature during PCR.')
    else:
        st.write(' The GC content of this sequence is within 40-65 range for ' +
                 'effective binding for PCR, which makes it a good candidate.')

stability(GCratio)

st.write("""

***
""")

#MRNA SEQUENCE
st.subheader('mRNA Sequence')
def transcription(DNAseq):
    for x in DNAseq:
       mRNA = DNAseq.replace('T', 'U')
    return mRNA

RNAseq = transcription(sequence)

st.text_area("DNA Sequence can be transcribed to the following mRNA sequence", RNAseq, height=400)
st.write("""

***
""")
# AMINO ACID SEQUENCE
st.subheader('Amino Acid Sequence')
def translation(RNA):
    codon_dictionary = {
                "AAA":"K", "AAC":"N", "AAG":"K", "AAU":"N",
                "ACA":"T", "ACC":"T", "ACG":"T", "ACU":"T",
                "AGA":"R", "AGC":"S", "AGG":"R", "AGU":"S",
                "AUA":"I", "AUC":"I", "AUG":"M", "AUU":"I",

                "CAA":"Q", "CAC":"H", "CAG":"Q", "CAU":"H",
                "CCA":"P", "CCC":"P", "CCG":"P", "CCU":"P",
                "CGA":"R", "CGC":"R", "CGG":"R", "CGU":"R",
                "CUA":"L", "CUC":"L", "CUG":"L", "CUU":"L",

                "GAA":"E", "GAC":"D", "GAG":"E", "GAU":"D",
                "GCA":"A", "GCC":"A", "GCG":"A", "GCU":"A",
                "GGA":"G", "GGC":"G", "GGG":"G", "GGU":"G",
                "GUA":"V", "GUC":"V", "GUG":"V", "GUU":"V",

                "UAA":"STOP", "UAC":"Y", "UAG":"STOP", "UAU":"T",
                "UCA":"S", "UCC":"S", "UCG":"S", "UCU":"S",
                "UGA":"STOP", "UGC":"C", "UGG":"W", "UGU":"C",
                "UUA":"L", "UUC":"F", "UUG":"L", "UUU":"F"
    }

    AAseq = ''

    for i in range(0, len(RNA) - (3 + len(RNA)%3), 3):
        if codon_dictionary[RNA[i:i+3]] == "STOP":
            break
        AAseq += codon_dictionary[RNA[i:i+3]]
    return AAseq

AAseq = translation(RNAseq)

st.text_area("mRNA Sequence can be transcribed to the following mRNA sequence", AAseq, height=400)

st.write("""

***
""")

#POLARITY
st.subheader('Polarity')
st.write('Each amino acid within the amino acid sequence carries specific polarity characteristics. Analysis ' +
         ' of this sequence, the amino acids are grouped by charge and polarity. This information can be useful' +
         ' in determining the interaction between the protein structure and its target.')

def polarity(amino_acid):

    polar_dict = {
            'A': 'nonpolar',
            'G': 'nonpolar',
            'V': 'nonpolar',
            'W': 'nonpolar',
            'I': 'nonpolar',
            'L': 'nonpolar',
            'M': 'nonpolar',
            'F': 'nonpolar',
            'P': 'nonpolar',
            'R': 'positive',
            'K': 'positive',
            'H': 'positive',
            'D': 'negative',
            'E': 'negative',
            'S': 'neutral',
            'T': 'neutral',
            'Y': 'neutral',
            'Q': 'neutral',
            'N': 'neutral',
            'C': 'neutral'
            }

    npcount = 0
    poscount = 0
    negcount = 0
    neutralcount = 0
    for o in amino_acid:
        if polar_dict[o] == "nonpolar":
            npcount += 1
            continue
        elif polar_dict[o] == "positive":
            poscount += 1
            continue

        elif polar_dict[o] == "negative":
            negcount += 1
            continue
        elif polar_dict[o] == "neutral":
            neutralcount += 1
            continue
    data = dict([("Nonpolar AA's", npcount), ("Polar-Positive AA's", poscount),
                 ("Polar-Negative AA's", negcount), ("Polar-Neutral AA's", neutralcount)])

    return data

data = polarity(AAseq)

# DATAFRAME

df1 = pd.DataFrame.from_dict(data, orient='index')
df1 = df1.rename({0: 'Count'}, axis='columns')
df1.reset_index(inplace=True)
df1 = df1.rename(columns={'index': 'Type'})
st.write(df1)

# POLARITY DISTRUBUTION BAR PLOT

st.write("""

    """)
lo = alt.Chart(df1).mark_bar().encode(
     x='Type',
     y='Count').configure_mark(opacity = 0.2, color = 'green')
lo = lo.properties(width=alt.Step(100))

st.write(lo)

st.write("""


***
""")

st.write("Special thanks to Dr. Katherine Herbert, Montclair University, for her continuous support during the development of this bioinformatics application. "
         + "There will be updates with new features in the new future. For any support please contact patels22@montclair.edu.")














