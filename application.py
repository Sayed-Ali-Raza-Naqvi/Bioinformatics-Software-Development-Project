import streamlit as st
import random
from Bio import pairwise2, SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqIO import parse, write
from io import StringIO
import matplotlib.pyplot as plt

# Function to calculate sequence length
def calculate_sequence_length(input_seq):
    if input_seq.startswith('>'):
        input_seq = input_seq.split('\n', 1)[1].replace('\n', '')
    return len(input_seq)

# Function to generate reverse complement
def generate_reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_complement = ''.join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement

# Function to align two sequences
def align_sequences(seq1, seq2, alignment_type):
    if alignment_type == "Global":
        alignments = pairwise2.align.globalxx(seq1, seq2)
    else:
        alignments = pairwise2.align.localxx(seq1, seq2)
    return alignments

# Function to randomize sequence
def randomize_sequence(seq):
    seq_list = list(seq)
    random.shuffle(seq_list)
    return ''.join(seq_list)

# Function to calculate GC content and create a plot
def calculate_gc_content(seq):
    gc_count = seq.count('G') + seq.count('C')
    gc_content = (gc_count / len(seq)) * 100
    return gc_content, gc_count, len(seq) - gc_count

# Function to find ORFs
def find_orfs(seq):
    orfs = []
    for strand, nuc in [(+1, seq), (-1, generate_reverse_complement(seq))]:
        for frame in range(3):
            trans = str(Seq(nuc[frame:]).translate(to_stop=True))
            if len(trans) >= 1:
                orfs.append(trans)
    return orfs

# Function to calculate isoelectric point
def calculate_isoelectric_point(protein_seq):
    analysed_seq = ProteinAnalysis(protein_seq)
    return analysed_seq.isoelectric_point()

# Function to predict protein secondary structure
def predict_secondary_structure(protein_seq):
    analysed_seq = ProteinAnalysis(protein_seq)
    return analysed_seq.secondary_structure_fraction()

# Function to convert file formats
def convert_file_format(input_file, input_format, output_format):
    input_handle = StringIO(input_file)
    output_handle = StringIO()
    sequences = list(parse(input_handle, input_format))

    if input_format == "fasta" and output_format == "fastq":
        for record in sequences:
            if "phred_quality" not in record.letter_annotations:
                record.letter_annotations["phred_quality"] = [40] * len(record.seq)
    write(sequences, output_handle, output_format)
    return output_handle.getvalue()

# Function to calculate genomic distance using Hamming distance
def calculate_hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        return "Sequences must be of the same length for Hamming distance"
    return sum(el1 != el2 for el1, el2 in zip(seq1, seq2))

# Function to calculate genomic distance using Jaccard distance
def calculate_jaccard_distance(seq1, seq2):
    set1 = set(seq1)
    set2 = set(seq2)
    return 1 - len(set1.intersection(set2)) / len(set1.union(set2))

# Home Page
def home():
    #st.title("Welcome to Enigma!")
    st.markdown("<p class='title-animation'></p>", unsafe_allow_html=True)
    st.markdown("<p style='font-size: 18px;'>A simple bioinformatics suite.</p>", unsafe_allow_html=True)
    

    # CSS animation for the title
    st.markdown("""
    <style>
    @keyframes typing-animation {
        0% { width: 0; }
        100% { width: 100%; }
    }

    .title-animation::after {
        content: 'Welcome to Enigma';
        display: inline-block;
        overflow: hidden;
        width: 0;
        animation: typing-animation 8s ease-in-out forwards;
        white-space: nowrap;
        font-size: 60px;
    }
    </style>
    """, unsafe_allow_html=True)

# Tools Page
def tools():
    st.title("Enigma Tools")
    st.write("Select a tool from below to use.")

    tool = st.radio("Select a tool", ["Sequence Length Calculator", "Sequence Reverse Complement Generator", "Sequence Alignment Viewer", "Sequence Randomizer", "GC Content Calculator", "ORF Finder", "Protein Isoelectric Point Calculator", "Protein Secondary Structure Predictor", "File Format Converter", "Genomic Distance Calculator"])

    if tool == "Sequence Length Calculator":
        # Sequence Length Calculator tool implementation
        st.header("Sequence Length Calculator")
        input_seq = st.text_area("Enter the sequence")
        if st.button("Calculate Sequence Length"):
            if input_seq:
                seq_length = calculate_sequence_length(input_seq)
                st.write(f"The length of the sequence is: {seq_length}")
            else:
                st.error("Please provide a sequence.")
    elif tool == "Sequence Reverse Complement Generator":
        # Sequence Reverse Complement Generator tool implementation
        st.header("Sequence Reverse Complement Generator")
        sequence = st.text_input("Enter a DNA sequence")
        if st.button("Generate Reverse Complement"):
            if sequence:
                reverse_complement = generate_reverse_complement(sequence)
                st.write(f"Reverse complement: {reverse_complement}")
            else:
                st.error("Please provide a sequence.")
    elif tool == "Sequence Alignment Viewer":
        # Sequence Alignment Viewer tool implementation
        st.header("Sequence Alignment Viewer")
        seq1 = st.text_area("Enter the first sequence", height=100)
        seq2 = st.text_area("Enter the second sequence", height=100)
        alignment_type = st.radio("Choose alignment type:", ("Global", "Local"))

        if st.button("Align Sequences"):
            if seq1 and seq2:
                alignments = align_sequences(seq1, seq2, alignment_type)
                st.write("Alignments:")
                for alignment in alignments:
                    st.text(pairwise2.format_alignment(*alignment))
            else:
                st.error("Please provide both sequences.")
    elif tool == "Sequence Randomizer":
        # Sequence Randomizer tool implementation
        st.header("Sequence Randomizer")
        random_seq = st.text_input("Enter a sequence to randomize")
        if st.button("Randomize Sequence"):
            if random_seq:
                randomized_sequence = randomize_sequence(random_seq)
                st.write(f"Randomized sequence: {randomized_sequence}")
            else:
                st.error("Please provide a sequence.")
    elif tool == "GC Content Calculator":
        # GC Content Calculator tool implementation
        st.header("GC Content Calculator")
        gc_seq = st.text_input("Enter a DNA sequence to calculate GC content")
        if st.button("Calculate GC Content"):
            if gc_seq:
                gc_content, gc_count, at_count = calculate_gc_content(gc_seq)
                st.write(f"GC content: {gc_content:.2f}%")
                
                fig, ax = plt.subplots()
                ax.pie([gc_count, at_count], labels=['GC', 'AT'], autopct='%1.1f%%', colors=['#66b3ff', '#ffcc99'])
                ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
                st.pyplot(fig)
            else:
                st.error("Please provide a sequence.")
    elif tool == "ORF Finder":
        # ORF Finder tool implementation
        st.header("ORF Finder")
        orf_seq = st.text_area("Enter a DNA sequence to find ORFs")
        if st.button("Find ORFs"):
            if orf_seq:
                orfs = find_orfs(orf_seq)
                st.write("Open Reading Frames (ORFs):")
                for i, orf in enumerate(orfs):
                    st.text(f"ORF {i+1}: {orf}")
            else:
                st.error("Please provide a sequence.")
    elif tool == "Protein Isoelectric Point Calculator":
        # Protein Isoelectric Point Calculator tool implementation
        st.header("Protein Isoelectric Point Calculator")
        protein_seq = st.text_input("Enter a protein sequence to calculate its isoelectric point")
        if st.button("Calculate Isoelectric Point"):
            if protein_seq:
                pi = calculate_isoelectric_point(protein_seq)
                st.write(f"Isoelectric point (pI): {pi:.2f}")
            else:
                st.error("Please provide a protein sequence.")
    elif tool == "Protein Secondary Structure Predictor":
        # Protein Secondary Structure Predictor tool implementation
        st.header("Protein Secondary Structure Predictor")
        secondary_structure_seq = st.text_input("Enter a protein sequence to predict its secondary structure")
        if st.button("Predict Secondary Structure"):
            if secondary_structure_seq:
                helix, sheet, coil = predict_secondary_structure(secondary_structure_seq)
                st.write(f"Helix: {helix*100:.2f}%")
                st.write(f"Sheet: {sheet*100:.2f}%")
                st.write(f"Coil: {coil*100:.2f}%")
            else:
                st.error("Please provide a protein sequence.")
    elif tool == "File Format Converter":
        # File Format Converter tool implementation
        st.header("File Format Converter")
        uploaded_file_format = st.file_uploader("Upload a file to convert", type=["fasta", "fastq"])
        input_format = st.selectbox("Select input format", ["fasta", "fastq"])
        output_format = st.selectbox("Select output format", ["fasta", "fastq"])
        if st.button("Convert File Format"):
            if uploaded_file_format:
                input_file_str = uploaded_file_format.read().decode("utf-8")
                converted_file_str = convert_file_format(input_file_str, input_format, output_format)
                st.download_button("Download converted file", data=converted_file_str, file_name=f"converted.{output_format}")
            else:
                st.error("Please upload a file.")
    elif tool == "Genomic Distance Calculator":
        # Genomic Distance Calculator tool implementation
        st.header("Genomic Distance Calculator")
        distance_seq1 = st.text_area("Enter the first sequence for distance calculation")
        distance_seq2 = st.text_area("Enter the second sequence for distance calculation")
        distance_method = st.selectbox("Select distance calculation method", ["Hamming", "Jaccard"])
        if st.button("Calculate Genomic Distance"):
            if distance_seq1 and distance_seq2:
                if distance_method == "Hamming":
                    distance = calculate_hamming_distance(distance_seq1, distance_seq2)
                    st.write(f"Hamming distance: {distance}")
                elif distance_method == "Jaccard":
                    distance = calculate_jaccard_distance(distance_seq1, distance_seq2)
                    st.write(f"Jaccard distance: {distance:.2f}")
            else:
                st.error("Please provide both sequences.")

# Pipeline Page
def pipeline():
    st.title("Enigma Pipeline")
    st.subheader("Sequence Input")

    sequence_input = ""
    input_method = st.radio("Choose input method:", ("Enter sequence", "Upload FASTA file"))

    if input_method == "Enter sequence":
        sequence_input = st.text_area("Enter DNA sequence")
    else:
        uploaded_file = st.file_uploader("Upload FASTA file", type=["fasta"])
        if uploaded_file is not None:
            sequence_data = uploaded_file.read().decode("utf-8")
            fasta_sequence = SeqIO.read(StringIO(sequence_data), "fasta")
            sequence_input = str(fasta_sequence.seq)

    if st.button("Submit"):
        if sequence_input:
            sequence = Seq(sequence_input.upper())

            # Transcription
            transcribed_seq = sequence.transcribe()

            # Translation
            translated_seq = transcribed_seq.translate()

            # Isoelectric Point Calculation
            protein_analysis = ProteinAnalysis(str(translated_seq))
            isoelectric_point = protein_analysis.isoelectric_point()

            # GC Content Calculation
            gc_content = calculate_gc_content(sequence)

            # ORF Finder
            orfs = find_orfs(str(sequence))

            # Display results
            st.subheader("Results")
            st.write("### Transcribed Sequence")
            st.write(transcribed_seq)

            st.write("### Translated Sequence")
            st.write(translated_seq)

            st.write("### Isoelectric Point of Protein")
            st.write(f"{isoelectric_point:.2f}")

            st.write("### GC Content")
            gc_content, gc_count, at_count = calculate_gc_content(sequence)
            st.write(f"GC content: {gc_content:.2f}%")
                
            fig, ax = plt.subplots()
            ax.pie([gc_count, at_count], labels=['GC', 'AT'], autopct='%1.1f%%', colors=['#66b3ff', '#ffcc99'])
            ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
            st.pyplot(fig)

            st.write("### Open Reading Frames (ORFs)")
            for i, orf in enumerate(orfs):
                st.text(f"ORF {i+1}: {orf}")
            
# About Page
def about():
    st.title("About Enigma")
    st.write("Enigma is a simple bioinformatics suite designed to provide various tools and utilities for bioinformatics analysis.")

    # Footer
    st.markdown("---")
    st.markdown("""
    <style>
    .footer {
        display: flex;
        position: fixed;
        bottom: 0;
        height: 70px; /* Adjust height as needed */
    }
    </style>
    """, unsafe_allow_html=True)
    st.markdown("<div class='footer'>Created by SAR Naqvi</div>", unsafe_allow_html=True)
    
# Streamlit app
st.set_page_config(page_title="Enigma", page_icon="🔬")

# Page Navigation
page = st.sidebar.radio("Go to", ["Home", "Tools", "Pipeline", "About"])

if page == "Home":
    home()
elif page == "Tools":
    tools()
elif page == "Pipeline":
    pipeline()
elif page == "About":
    about()
           
