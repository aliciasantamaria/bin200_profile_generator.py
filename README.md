# generate_200bp_binary_profiles

This Python script processes genomic peak data (in BED format) and converts binary signal files (e.g. from ChIP-seq or similar experiments) into **200 bp window-based binary genomic profiles**.  
It can handle both **female (F)** and **male (M)** datasets, generating one binary output file per sample.

---

## 🧬 Overview

The script:
1. Reads a BED file containing genomic peaks.
2. Divides the chromosome into **200 bp non-overlapping windows**.
3. Reads binary signal files for each sample (e.g. `*_chr10_binary.txt`).
4. Maps each binary value (0/1) from peaks to the corresponding 200 bp window.
5. Generates binary output files, preserving the input file headers.

This approach enables downstream analyses of epigenetic or chromatin accessibility patterns at fixed resolution (200 bp), facilitating comparison across samples or groups.

---

## ⚙️ Requirements

- Python ≥ 3.8  
- [pandas](https://pandas.pydata.org/)  
- [pybedtools](https://daler.github.io/pybedtools/)  
- [tqdm](https://tqdm.github.io/)  

You can install all dependencies via:

```bash
pip install pandas pybedtools tqdm
```

---

## 🚀 Usage

```bash
python generate_200bp_binary_profiles.py peaks.bed \
  --input_dir_F your/input/data/female \
  --input_dir_M your/input/data/male \
  --output_dir_F your/output/folder/female \
  --output_dir_M your/output/folder/male
```

### **Arguments**

| Argument | Description |
|-----------|-------------|
| `peaks_file` | Input BED file containing peaks (e.g. `chr10_peaks.bed`). |
| `--input_dir_F` | Directory with binary input files for female samples. |
| `--input_dir_M` | Directory with binary input files for male samples. |
| `--output_dir_F` | Output directory for processed female files. |
| `--output_dir_M` | Output directory for processed male files. |

All directories will be created automatically if they don’t exist.

---

## 📂 Input file format

Each input binary file should follow this pattern:

```
# Sample information (header line 1)
# Metadata (header line 2)
1
0
0
1
...
```

- The filename must include the chromosome (e.g. `_chr10_binary.txt`).
- The number of binary values must match the number of peaks in the BED file.

---

## 💾 Output

For each sample, the script generates a new binary file with the same name, but containing binary values corresponding to **200 bp windows** instead of individual peaks.

Output structure:

```
# Sample information
# Metadata
0
1
0
...
```

Each line corresponds to one 200 bp genomic window.

---

## 🧠 Example

If your input BED (`chr10_peaks.bed`) defines peaks along chromosome 10, and your input directories contain files like:

```
outputF/sample1_chr10_binary.txt
outputF/sample2_chr10_binary.txt
outputM/sample3_chr10_binary.txt
```

Then the script will produce:

```
output200F_allchr/sample1_chr10_binary.txt
output200F_allchr/sample2_chr10_binary.txt
output200M_allchr/sample3_chr10_binary.txt
```

Each output file will contain the binary signal per 200 bp window.

---

## ⚖️ License

This project is released under the **MIT License** — see the [LICENSE](LICENSE) file for details.

```
MIT License

Copyright (c) 2025 [Your Name]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
```

---

## ✨ Citation / Credits

If you use this script in your analysis or publication, please cite it as:

> Santamaría-Quiles, Alicia, *generate_200bp_binary_profiles.py*, GitHub repository (2025)

---

## 🧩 Author

**Alicia Santamaría Quiles**  
📧 [alicia.santamaria@cabimer.es]


