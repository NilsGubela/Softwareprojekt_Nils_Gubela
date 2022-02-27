# Pausing site finding algorithm

The present program is a Python wrapper for `DrTransformer` that given a RNA sequence structure pair searches for pausing sites to increase the structures occupancy at the end of transcription.
`DrTransformer` and Python 3 must be installed separately. See below for an example workflow.

## Purpose

The pausing site finding algorithm takes a list of sequence structure pairs as input and searches for pausing sites to increase the structures occupancy at the end of transcription. The output is a suggestion of pausing sites in the format "site=length" that can be used as input for DrTransformer
```sh
    echo ACUUAGUGGACUUCCAUAAACAAAGUCUCCACAAUAUUGGCACUUCCCGG | DrTransformer --pause-site 29=5
```
The generated csv file contains additional information such as the initial occupancy of the strucutre, the increased occupancy of the structure and in case the structure is not observable in the DrTransformer output it will list the closest observable structure w.r.t. base pair distance.

Please contact the [developer documentation](https://github.com/NilsGubela/Softwareprojekt_Nils_Gubela/blob/main/docs/dev_doc.pdf) of the [project report](https://github.com/NilsGubela/Softwareprojekt_Nils_Gubela/blob/main/docs/Report.pdf) for further information.

## Flags

 * -h, --help: Shows help message
 * -m, --mode: Either "adpative" or "exhaustive". The string "adaptive" stands for adaptive walk and "exhaustive" for exhaustive search and refers to the method of search that is used to find pausing sites. Default is adaptive walk.
 * -i, --input: Name of input file. Input file needs to be text file. Suffix ".txt" may not be hand over.
 * -o, --output: Name of output file. Comma seperated output file will be written. Suffix ".csv" may not be hand over. Default output file name is "output".
 * --o-prune: Flag for DrTransformer call. Will be used for any call. "Occupancy threshold to prune structures from the network. The structures with lowest occupancy are
                        removed until at most o-prune occupancy has beem removed from the total population. (default: 0.1)"
 * -t, --temp: Flag for DrTransformer call. Will be used for any call. Rescale energy parameters to a temperature of temp C. (default: 37.0). 

## Example

 * Create a text file that contains a sequence structure pair. There needs to be a space character between the sequence and the structure:
   ```sh
    echo "ACUUAGUGGACUUCCAUAAACAAAGUCUCCACAAUAUUGGCACUUCCCGG .....((((....)))).....((((..(((......))).))))....." > input.txt
    ```
 * Call the pausing site finding algorithm, specify the input file, the mode of search and the output 
file:
   ```sh
   python pause_finding.py -i input.txt -o output -m exhaustive 
    ```
 * Find the result as print out in the terminal and more detailed information in the output file.
 * alternatively the pruning parameter of DrTransformer can be set manually or the energy parameters can be rescaled to a different temperature:
 ```sh
   python pause_finding.py -i input.txt -o output_prune -m exhaustive --o-prune 0.001
   python pause_finding.py -i input.txt -o output_temp -m exhaustive -t 10
 ```
 

## Installation
Please install Python 3 and `DrTransformer` [v0.9](https://github.com/bad-ants-fleet/ribolands/tree/development). Navigate to the directory that contains the file pause_finding.py.

