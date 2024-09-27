import tkinter as tk
from tkinter import filedialog, messagebox
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


class GCVirolensApp:
    def __init__(self, root):
        self.root = root
        self.root.title("GCVirolens")

        # Label for the application name
        self.label = tk.Label(root, text="GCVirolens", font=("Arial", 24),fg="blue")
        self.label.grid(column=0, row=0, columnspan=4, pady=10)

        self.label2 = tk.Label(root, text="Gene-wise GC content\nanalysis of Virus Genomes", font=("Arial", 15))
        self.label2.grid(column=0, row=1, columnspan=4, pady=10)

        # Button to upload a FASTA file
        self.fasta_button = tk.Button(root, text="Upload genome file i.e  .fasta format",relief="sunken", command=self.upload_fasta, borderwidth=3,bg= "bisque")
        self.fasta_button.grid(column=0, row=2, columnspan=2, pady=10, padx=20)

        # Button to upload a GFF file
        self.gff_button = tk.Button(root, text="Upload annotation file i.e  .gff format",relief="sunken", command=self.upload_gff, borderwidth=3,bg= "bisque")
        self.gff_button.grid(column=2, row=2, columnspan=2, pady=10, padx=20)

        self.fasta_label = tk.Label(root, text="No file selected", font=("Arial", 10))
        self.fasta_label.grid(column=0, row=3, columnspan=2, pady=10)

        self.gff_label = tk.Label(root, text="No file selected", font=("Arial", 10))
        self.gff_label.grid(column=2, row=3, columnspan=2, pady=10)

        # Button to execute GC content analysis
        self.analyze_button = tk.Button(root, text="Execute GC Content Analysis",relief="raised", command=self.execute_analysis, borderwidth=3,bg="orange")
        self.analyze_button.grid(column=1, row=4, columnspan=2, pady=10)

        self.exec_label = tk.Label(root, text="Analysis yet to be initiated", font=("Arial", 10))
        self.exec_label.grid(column=1, row=5, columnspan=2, pady=10)

        # Button to save the result file
        self.save_button = tk.Button(root, text="Save Result",relief="raised", command=self.save_result, borderwidth=3,bg="turquoise")
        self.save_button.grid(column=1, row=6, columnspan=2, pady=10)

        self.result_label = tk.Label(root, text="No results to save", font=("Arial", 10))
        self.result_label.grid(column=1, row=7, columnspan=2, pady=10)

        self.contri_label = tk.Label(root, text="Developement Team:\n\nSumeet Kumar Parida,\nCenter for Post-Graduate Studies,\nOUAT, Bhubaneswar, Odisha,\nIndia", font=("Arial", 8), justify='left')
        self.contri_label.grid(column=0, row=8, columnspan=2, pady=10, sticky='w', padx=20)
        
        self.contri_label = tk.Label(root, text="Dr. Samarth Godara, \nICAR-IASRI,\nNew Delhi,\nIndia", font=("Arial", 8), justify='left')
        self.contri_label.grid(column=1, row=8, columnspan=2, pady=10, sticky='s', padx=20)
        
        self.contri_label = tk.Label(root, text="\n\nDr. Shbana Begam,\nICAR-NIPB,\nNew Delhi,\nIndia", font=("Arial", 8), justify='left')
        self.contri_label.grid(column=2, row=8, columnspan=2, pady=10, sticky='e', padx=20)
        


        # Variables to store file paths
        self.fasta_file = None
        self.gff_file = None
        self.result = None

    def upload_fasta(self):
        self.fasta_file = filedialog.askopenfilename(title="Select FASTA File", filetypes=(("FASTA files", "*.fasta *.fa *.fna"), ("All files", "*.*")))
        if self.fasta_file:
            self.fasta_label.config(text=self.fasta_file.split("/")[-1])
            messagebox.showinfo("File Selected", f"FASTA file selected: {self.fasta_file}")
        else:
            messagebox.showwarning("No File", "No FASTA file selected")

    def upload_gff(self):
        self.gff_file = filedialog.askopenfilename(title="Select GFF File", filetypes=(("GFF files", "*.gff"), ("All files", "*.*")))
        if self.gff_file:
            self.gff_label.config(text=self.gff_file.split("/")[-1])
            messagebox.showinfo("File Selected", f"GFF file selected: {self.gff_file}")
        else:
            messagebox.showwarning("No File", "No GFF file selected")

    def execute_analysis(self):
        if not self.fasta_file or not self.gff_file:
            messagebox.showerror("Missing Files", "Please upload both FASTA and GFF files before analysis.")
            return

        try:
            recs = []
            with open(self.gff_file, "r") as f:
                for line in f:
                    if line[0] != '#':
                        att = line.strip().split('\t')
                        recs.append(att)

            chr_list = []
            type_list = []
            start_list = []
            end_list = []
            id_list = []

            # Extract the relevant information from the records
            for rec in recs:
                chr_list.append(rec[0])
                type_list.append(rec[2])
                start_list.append(int(rec[3]) - 1)  # Convert to 0-based indexing
                end_list.append(int(rec[4]))

                # Extracting the gene ID from the Attributes column
                attributes = rec[8]
                gene_id = ''
                for attribute in attributes.split(';'):
                    if attribute.strip().startswith('ID='):
                        gene_id = attribute.strip().split('=')[1]
                        break
                id_list.append(gene_id)

            gff_data = pd.DataFrame({
                'genbank_id': chr_list,
                'type': type_list,
                'start': start_list,
                'end': end_list,
                'g_id': id_list
            })

            # Extract gene records
            gene_df = gff_data[gff_data["type"] == "gene"]

            # Load the sequences from the FASTA file
            sequences = []
            for record in SeqIO.parse(self.fasta_file, "fasta"):
                sequences.append(str(record.seq))

            # Assuming you want to work with the first sequence from the FASTA file
            sequence = sequences[0]

            gc_list = []

            # Calculate GC content for each gene
            for i in range(len(gene_df)):
                start_pos = gene_df["start"].values[i]
                end_pos = gene_df["end"].values[i]
                gene_seq = sequence[start_pos:end_pos]
                gc_content = gc_fraction(gene_seq) * 100
                gc_list.append(gc_content)

            # Prepare the output DataFrame
            output = pd.DataFrame({
                "g_id": gene_df["g_id"],
                "type": gene_df["type"],
                "start": gene_df["start"],
                "end": gene_df["end"],
                "GC_content": gc_list
            })

            # Store the result for later saving
            self.result = output.to_csv(index=False)
            self.exec_label.config(text="GC Content Analysis Complete")
            messagebox.showinfo("Analysis Complete", "GC content analysis is complete.")

        except Exception as e:
            messagebox.showerror("Error", f"An error occurred during analysis: {e}")

    def save_result(self):
        if self.result:
            save_path = filedialog.asksaveasfilename(title="Save Result", defaultextension=".csv", filetypes=(("CSV files", "*.csv"), ("All files", "*.*")))
            if save_path:
                with open(save_path, 'w') as file:
                    file.write(self.result)
                self.result_label.config(text=f"Result saved to: {save_path.split('/')[-1]}")
                messagebox.showinfo("File Saved", f"Result saved to: {save_path}")
            else:
                messagebox.showwarning("No File", "Save operation canceled")
        else:
            messagebox.showerror("No Result", "No analysis result to save. Please execute the analysis first.") 


if __name__ == "__main__":
    root = tk.Tk()
    app = GCVirolensApp(root)
    
    root.grid_rowconfigure(2, weight=0)
    root.grid_rowconfigure(4, weight=0)
    root.grid_rowconfigure(6, weight=0)
    root.grid_rowconfigure(8, weight=0)
    root.grid_columnconfigure(0, weight=0)
    root.grid_columnconfigure(1, weight=0)
    root.grid_columnconfigure(2, weight=0)
    
    root.mainloop()