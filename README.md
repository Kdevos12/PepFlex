# **PepFlex: A Modular Framework for Peptide Screening and Evolution**

PepFlex is a Python module designed to provide a flexible and extensible framework for in silico peptide screening and evolutionary studies. It offers tools for managing amino acid data, assembling peptides, applying various genetic operations (mutations and crossovers), and processing evolutionary rounds through a customizable evaluation pipeline.  
**Note on RDKit Dependency:** PepFlex utilizes RDKit for chemical functionalities, including SMILES canonicalization and peptide assembly. If RDKit is not installed, these functionalities will be limited, and a warning will be issued upon module import.

### **1\. Amino Acid DataFrame Utilities and Conversions**

PepFlex maintains a global DataFrame, AMINO\_ACID\_DF, which stores information about standard amino acids, including their one-letter, three-letter, and SMILES representations. This DataFrame can be extended to include non-canonical or modified amino acids.

* **AMINO\_ACID\_DF (Global DataFrame)**  
  * **Description:** A pandas DataFrame indexed by the 3-letter code, containing 1L (one-letter code), SMILES (original SMILES string), and SMILES\_Canon (canonical RDKit SMILES).  
  * **Purpose:** Serves as the central repository for amino acid data, facilitating conversions and lookups.  
* **\_canonicalize\_smiles(smiles: str) \-\> str**  
  * **Description:** Converts a SMILES string to its canonical RDKit form for reliable comparison. This is an internal helper function.  
  * **Arguments:**  
    * smiles (str): The SMILES string to canonicalize.  
  * **Returns:** str: The canonical SMILES string, or an empty string if invalid or RDKit is not available.  
* **add\_amino\_acid(three\_letter\_code: str, one\_letter\_code: str, smiles: str)**  
  * **Description:** Adds a new amino acid or updates an existing one in the AMINO\_ACINO\_DF using its 3-letter code as the primary key.  
  * **Arguments:**  
    * three\_letter\_code (str): The 3-letter code (e.g., "Ala", "Cys").  
    * one\_letter\_code (str): The 1-letter code (e.g., "A", "C").  
    * smiles (str): The SMILES string representation.  
  * **Behavior:** Prints a warning if the code already exists and updates it. Prints an error if the provided SMILES is invalid.  
* **get\_smiles\_from\_3L(three\_letter\_code: str) \-\> str**  
  * **Description:** Retrieves the original SMILES for a given 3-letter code.  
  * **Arguments:**  
    * three\_letter\_code (str): The 3-letter code of the amino acid.  
  * **Returns:** str: The SMILES string.  
  * **Raises:** ValueError: If the 3-letter code is not found.  
* **get\_1L\_from\_3L(three\_letter\_code: str) \-\> str**  
  * **Description:** Retrieves the 1-letter code for a given 3-letter code.  
  * **Arguments:**  
    * three\_letter\_code (str): The 3-letter code of the amino acid.  
  * **Returns:** str: The 1-letter code.  
  * **Raises:** ValueError: If the 3-letter code is not found.  
* **get\_smiles\_from\_1L(one\_letter\_code: str) \-\> str**  
  * **Description:** Retrieves the SMILES for a given 1-letter code.  
  * **Arguments:**  
    * one\_letter\_code (str): The 1-letter code of the amino acid.  
  * **Returns:** str: The SMILES string.  
  * **Raises:** ValueError: If the 1-letter code is not found.  
* **get\_3L\_from\_1L(one\_letter\_code: str) \-\> str**  
  * **Description:** Retrieves the 3-letter code for a given 1-letter code.  
  * **Arguments:**  
    * one\_letter\_code (str): The 1-letter code of the amino acid.  
  * **Returns:** str: The 3-letter code.  
  * **Raises:** ValueError: If the 1-letter code is not found.  
* **get\_1L\_from\_smiles(smiles: str) \-\> str**  
  * **Description:** Retrieves the 1-letter code for a given SMILES, using its canonical form for lookup.  
  * **Arguments:**  
    * smiles (str): The SMILES string of the amino acid.  
  * **Returns:** str: The 1-letter code.  
  * **Raises:** ValueError: If the SMILES is invalid or not found.  
* **get\_3L\_from\_smiles(smiles: str) \-\> str**  
  * **Description:** Retrieves the 3-letter code for a given SMILES, using its canonical form for lookup.  
  * **Arguments:**  
    * smiles (str): The SMILES string of the amino acid.  
  * **Returns:** str: The 3-letter code.  
  * **Raises:** ValueError: If the SMILES is invalid or not found.  
* **save\_amino\_acid\_data(file\_path: str, file\_format: str \= 'csv')**  
  * **Description:** Saves the AMINO\_ACID\_DF to a specified file.  
  * **Arguments:**  
    * file\_path (str): The path to save the file.  
    * file\_format (str): The format to save the file ('csv' or 'parquet'). Defaults to 'csv'.  
* **load\_amino\_acid\_data(file\_path: str, file\_format: str \= 'csv', overwrite: bool \= True)**  
  * **Description:** Loads amino acid data from a file and updates the global AMINO\_ACID\_DF.  
  * **Arguments:**  
    * file\_path (str): The path to the file.  
    * file\_format (str): The format of the file ('csv' or 'parquet'). Defaults to 'csv'.  
    * overwrite (bool): If True, overwrites the existing AMINO\_ACID\_DF. If False, merges the loaded data, keeping the last duplicate entry.

### **2\. Peptide Assembler (Beta)**

This function allows the construction of a complete peptide molecule from a sequence of amino acid SMILES strings.

* **assemble\_peptide(residue\_list: Union\[List\[str\], str\], input\_type: int \= 0\) \-\> Optional\[Chem.Mol\]**  
  * **Description:** Assembles a peptide from a list of amino acid SMILES strings (N- to C-terminus order) or a hyphen-separated string of 3-letter codes. Each input SMILES must represent an amino acid with free N and C termini (e.g., N\[C@@H\](R)C(=O)O).  
  * **Arguments:**  
    * residue\_list (Union\[List\[str\], str\]): A list of SMILES strings or a string of 3-letter codes (e.g., "Ala-Glu-Asp-Gly").  
    * input\_type (int): 0 for a list of SMILES (default), 1 for a 3-letter code string.  
  * **Returns:** Optional\[Chem.Mol\]: An RDKit Mol object representing the assembled peptide, or None if the input list is empty.  
  * **Raises:**  
    * ValueError: If RDKit is not installed, SMILES cannot be parsed, or a 3-letter code is not found.  
    * RuntimeError: If peptide bond formation fails between residues.

### **3\. Miscellaneous Utilities**

This section includes functions for converting peptide sequences between different amino acid representations and a function for performing genetic crossover operations.

* **convert\_1L\_to\_smiles\_list(one\_letter\_str: str) \-\> List\[str\]**  
  * **Description:** Converts a string of 1-letter amino acid codes to a list of their corresponding SMILES strings.  
  * **Arguments:**  
    * one\_letter\_str (str): A string of 1-letter codes (e.g., "AEDG").  
  * **Returns:** List\[str\]: A list of SMILES strings.  
* **convert\_smiles\_list\_to\_3L(smiles\_list: List\[str\]) \-\> str**  
  * **Description:** Converts a list of SMILES strings to a hyphen-separated string of 3-letter amino acid codes.  
  * **Arguments:**  
    * smiles\_list (List\[str\]): A list of SMILES strings.  
  * **Returns:** str: A string of 3-letter codes (e.g., "Ala-Glu-Asp-Gly").  
* **convert\_1L\_str\_to\_3L(one\_letter\_str: str) \-\> str**  
  * **Description:** Converts a 1-letter amino acid string to a hyphen-separated 3-letter amino acid string. This only works for canonical amino acids present in the global DataFrame.  
  * **Arguments:**  
    * one\_letter\_str (str): A string of 1-letter codes (e.g., "AEDG").  
  * **Returns:** str: A string of 3-letter codes (e.g., "Ala-Glu-Asp-Gly").  
* **convert\_3L\_str\_to\_1L(three\_letter\_str: str) \-\> str**  
  * **Description:** Converts a hyphen-separated 3-letter amino acid string to a 1-letter amino acid string. This only works for canonical amino acids present in the global DataFrame.  
  * **Arguments:**  
    * three\_letter\_str (str): A string of 3-letter codes (e.g., "Ala-Glu-Asp-Gly").  
  * **Returns:** str: A string of 1-letter codes (e.g., "AEDG").  
* **convert\_smiles\_list\_to\_1L(smiles\_list: List\[str\]) \-\> str**  
  * **Description:** Converts a list of SMILES strings to a 1-letter amino acid string. This only works for canonical amino acids present in the global DataFrame.  
  * **Arguments:**  
    * smiles\_list (List\[str\]): A list of SMILES strings.  
  * **Returns:** str: A string of 1-letter amino acid codes (e.g., "AEDG").  
* **peptide\_crossover(parent1\_smiles: List\[str\], parent2\_smiles: List\[str\]) \-\> Tuple\[List\[str\], List\[str\]\]**  
  * **Description:** Performs a single-point crossover operation between two parent peptides. The cut-off point is determined as the middle of each peptide, with a random choice for odd-length peptides. The segments are then exchanged to form two child peptides. Handles empty or very short peptides.  
  * **Arguments:**  
    * parent1\_smiles (List\[str\]): The SMILES sequence of the first parent.  
    * parent2\_smiles (List\[str\]): The SMILES sequence of the second parent.  
  * **Returns:** Tuple\[List\[str\], List\[str\]\]: A tuple containing the SMILES sequences of the two child peptides.

### **4\. Peptide Mutator**

The PeptideMutator class enables the application of various mutation rules to peptide sequences. It supports predefined mutation types and allows for custom user-defined mutation functions.

* **class PeptideMutator**  
  * **\_\_init\_\_(self)**  
    * **Description:** Initializes the PeptideMutator with an empty list of mutation rules and sets the default amino acid DataFrame.  
  * **add\_mutation\_rule(self, mutation\_type: Union\[str, Callable\], probability: float, \*\*kwargs)**  
    * **Description:** Adds a mutation rule to the mutator.  
    * **Arguments:**  
      * mutation\_type (Union\[str, Callable\]): The type of mutation. Can be a string for predefined rules or a callable function for custom rules.  
        * **Predefined Types:** 'intra\_sequence\_suppression', 'n\_terminal\_suppression', 'c\_terminal\_suppression', 'n\_terminal\_addition', 'c\_terminal\_addition', 'intra\_sequence\_addition', 'inter\_mutation', 'n\_terminal\_mutation', 'c\_terminal\_mutation', 'segment\_duplication', 'segment\_inversion', 'segment\_translocation'.  
        * **Custom Functions:** Should accept smiles\_list: List\[str\] as the first argument and \*\*kwargs for additional parameters, returning a List\[str\].  
      * probability (float): The probability (0.0 to 1.0) of applying this rule.  
      * \*\*kwargs: Additional parameters specific to the mutation rule. For predefined rules, these might include allowed\_amino\_acids\_df, source\_amino\_acids\_df, target\_amino\_acids\_df, min\_length, max\_length.  
    * **Raises:**  
      * ValueError: If probability is out of range or mutation type is unknown.  
      * TypeError: If DataFrame parameters are not pandas.DataFrame.  
  * **apply\_mutations(self, initial\_smiles\_list: list) \-\> list**  
    * **Description:** Applies all configured mutation rules to a given peptide SMILES sequence based on their defined probabilities.  
    * **Arguments:**  
      * initial\_smiles\_list (List\[str\]): The initial SMILES sequence of the peptide.  
    * **Returns:** List\[str\]: The mutated SMILES sequence.

### **5\. Peptide Class**

The Peptide class encapsulates all information related to a single peptide, including its sequence, unique identifier, creation date, and a history of events.

* **class Peptide**  
  * **\_\_init\_\_(self, smiles\_sequence: list, peptide\_id: str \= None, creation\_date: datetime \= None, history: list \= None, source\_generation\_params: dict \= None)**  
    * **Description:** Initializes a Peptide object.  
    * **Arguments:**  
      * smiles\_sequence (list): A list of SMILES strings representing the amino acid sequence.  
      * peptide\_id (str, optional): Unique identifier for the peptide. Generated if None.  
      * creation\_date (datetime, optional): Timestamp of peptide creation. Defaults to now.  
      * history (list, optional): A list of dictionaries detailing events in the peptide's lifecycle.  
      * source\_generation\_params (dict, optional): Parameters related to how the peptide was generated.  
  * **set\_sequences\_from\_global\_converter(self)**  
    * **Description:** Sets the 1-letter and 3-letter sequence representations using global conversion functions. Handles cases where SMILES might not be recognized by marking them as "N/A" and logging.  
  * **one\_letter\_sequence (property)**  
    * **Description:** Returns the 1-letter code sequence of the peptide.  
  * **three\_letter\_sequence (property)**  
    * **Description:** Returns the 3-letter code sequence of the peptide.  
  * **add\_history\_entry(self, event\_type: str, details: dict)**  
    * **Description:** Adds an entry to the peptide's history log.  
    * **Arguments:**  
      * event\_type (str): The type of event (e.g., "mutation", "crossover", "replenishment").  
      * details (dict): A dictionary containing event-specific details.  
  * **to\_dict(self) \-\> Dict**  
    * **Description:** Converts the Peptide object to a dictionary representation.  
    * **Returns:** Dict: A dictionary containing all peptide attributes.  
  * **from\_dict(cls, data: Dict) \-\> 'Peptide' (classmethod)**  
    * **Description:** Creates a Peptide object from a dictionary representation.  
    * **Arguments:**  
      * data (Dict): A dictionary containing peptide data.  
    * **Returns:** Peptide: A new Peptide instance.  
  * **\_\_repr\_\_(self) \-\> str**  
    * **Description:** Provides a string representation of the Peptide object for debugging.  
  * **\_\_eq\_\_(self, other: object) \-\> bool**  
    * **Description:** Defines equality comparison based on peptide\_id.  
  * **\_\_hash\_\_(self) \-\> int**  
    * **Description:** Defines the hash for the object based on peptide\_id.

### **6\. PeptidePoolManager Class**

The PeptidePoolManager class is responsible for managing a collection of Peptide objects, providing functionalities to add, retrieve, remove, and persist the peptide pool.

* **class PeptidePoolManager**  
  * **\_\_init\_\_(self)**  
    * **Description:** Initializes an empty dictionary to store Peptide objects, keyed by their peptide\_id.  
  * **add\_peptide(self, peptide: Peptide)**  
    * **Description:** Adds a Peptide object to the pool.  
    * **Arguments:**  
      * peptide (Peptide): The Peptide object to add.  
    * **Raises:** TypeError: If the object is not a Peptide instance.  
  * **get\_peptide(self, peptide\_id: str) \-\> Optional\[Peptide\]**  
    * **Description:** Retrieves a peptide by its ID.  
    * **Arguments:**  
      * peptide\_id (str): The ID of the peptide.  
    * **Returns:** Optional\[Peptide\]: The Peptide object if found, else None.  
  * **get\_all\_peptides(self) \-\> List\[Peptide\]**  
    * **Description:** Returns a list of all Peptide objects currently in the pool.  
    * **Returns:** List\[Peptide\]: A list of Peptide objects.  
  * **remove\_peptide(self, peptide\_id: str)**  
    * **Description:** Removes a peptide from the pool by its ID.  
    * **Arguments:**  
      * peptide\_id (str): The ID of the peptide to remove.  
  * **clear\_pool(self)**  
    * **Description:** Clears all peptides from the pool.  
  * **get\_pool\_size(self) \-\> int**  
    * **Description:** Returns the current number of peptides in the pool.  
    * **Returns:** int: The number of peptides.  
  * **to\_dataframe(self) \-\> pd.DataFrame**  
    * **Description:** Converts the peptide pool into a pandas DataFrame. History and source generation parameters are serialized to JSON strings for storage.  
    * **Returns:** pd.DataFrame: A DataFrame representation of the pool.  
  * **from\_dataframe(self, df: pd.DataFrame)**  
    * **Description:** Loads peptides into the pool from a pandas DataFrame. History and source generation parameters are deserialized from JSON strings.  
    * **Arguments:**  
      * df (pd.DataFrame): The DataFrame containing peptide data.  
  * **save\_pool(self, file\_path: str, file\_format: str \= 'parquet')**  
    * **Description:** Saves the current peptide pool to a file.  
    * **Arguments:**  
      * file\_path (str): The path to save the file.  
      * file\_format (str): The format to save the file ('parquet', 'csv', or 'json'). Defaults to 'parquet'.  
  * **load\_pool(self, file\_path: str, file\_format: str \= 'parquet')**  
    * **Description:** Loads a peptide pool from a file into the current manager.  
    * **Arguments:**  
      * file\_path (str): The path to the file to load.  
      * file\_format (str): The format of the file ('parquet', 'csv', or 'json'). Defaults to 'parquet'.

### **7\. Evaluator**

The Evaluator class manages the evaluation pipeline for peptides. It converts peptides to DataFrames, applies user-defined processing steps (En), and ranks them using a specified ranker function. Intermediate DataFrames can be saved for logging and debugging.

* **class Evaluator**  
  * **\_\_init\_\_(self, evaluation\_pipeline: List\[Callable\[\[pd.DataFrame\], pd.DataFrame\]\], ranker\_function: Callable\[\[pd.DataFrame\], Tuple\[List\['Peptide'\], Optional\[List\[float\]\]\]\]):**  
    * **Description:** Initializes the Evaluator with a pipeline of evaluation steps and a ranker function.  
    * **Arguments:**  
      * evaluation\_pipeline (List\[Callable\[\[pd.DataFrame\], pd.DataFrame\]\]): A list of functions (En steps). Each function takes a DataFrame and returns a DataFrame. These functions define the sequence of data processing and feature generation.  
      * ranker\_function (Callable\[\[pd.DataFrame\], Tuple\[List\['Peptide'\], Optional\[List\[float\]\]\]\]): A function that takes the final processed DataFrame and returns a tuple: (List\[Peptide\] of ranked peptides, Optional\[List\[float\]\] of their corresponding scores). This function is responsible for the final selection/ordering of peptides.  
  * **\_e0\_to\_dataframe(self, peptides: List\['Peptide'\]) \-\> pd.DataFrame**  
    * **Description:** Converts a list of Peptide objects into the initial DataFrame (E0 step). Each row represents a peptide, with its full data. This is an internal helper function.  
    * **Arguments:**  
      * peptides (List\['Peptide'\]): The list of Peptide objects.  
    * **Returns:** pd.DataFrame: The initial DataFrame for evaluation.  
  * **\_save\_dataframe(self, df: pd.DataFrame, step\_name: str, round\_identifier: str \= None)**  
    * **Description:** Saves an intermediate DataFrame to a CSV file for logging and debugging. Filename includes peptide count, feature count, round identifier, step name, and a hash. This is an internal helper function.  
    * **Arguments:**  
      * df (pd.DataFrame): The DataFrame to save.  
      * step\_name (str): The name of the evaluation step (for filename).  
      * round\_identifier (str, optional): An identifier for the current round (for filename).  
  * **evaluate(self, peptides: List\['Peptide'\], round\_identifier: str \= None) \-\> Tuple\[List\['Peptide'\], Optional\[List\[float\]\]\]**  
    * **Description:** Executes the evaluation pipeline on a list of peptides and returns the ranked peptides.  
    * **Arguments:**  
      * peptides (List\['Peptide'\]): The list of Peptide objects to evaluate.  
      * round\_identifier (str, optional): An identifier for the current round, used for logging.  
    * **Returns:** Tuple\[List\['Peptide'\], Optional\[List\[float\]\]\]: A tuple containing:  
      * A list of Peptide objects, ranked according to the ranker\_function.  
      * An optional list of scores corresponding to the ranked peptides (used for truncation).

### **8\. Round Processor (Evolution Logic)**

The PoolRoundProcessor orchestrates a single round of peptide evolution, applying a sequence of mutation, crossover, evaluation, replenishment, and truncation steps as defined by its pipeline.

* **class PoolRoundProcessor**  
  * **\_\_init\_\_(self)**  
    * **Description:** Initializes the round processor with an empty pipeline and no generation function set.  
  * **add\_pipeline\_step(self, step\_type: str, step\_function: Callable, name: str, \*\*params)**  
    * **Description:** Adds a step to the round processing pipeline.  
    * **Arguments:**  
      * step\_type (str): The type of step ('mutation', 'crossover', 'evaluation', 'replenishment', 'truncation', 'custom').  
      * step\_function (Callable): The function to execute for this step. Its signature depends on step\_type.  
      * name (str): A descriptive name for the step.  
      * \*\*params: Parameters specific to the step\_function.  
    * **Raises:** TypeError: If step\_function is not callable.  
  * **set\_generation\_function(self, gen\_func: Callable\[\[int\], list\['Peptide'\]\])**  
    * **Description:** Sets a function to generate new peptides if the pool is too small (used by the replenishment step).  
    * **Arguments:**  
      * gen\_func (Callable\[\[int\], list\['Peptide'\]\]): A function that takes the number of peptides to generate and returns a list of Peptide objects.  
  * **\_execute\_mutation\_step(self, peptides: List\['Peptide'\], mutator: 'PeptideMutator', probability\_of\_application: float) \-\> List\['Peptide'\]**  
    * **Description:** Internal helper to apply mutations to a list of peptides.  
  * **\_execute\_crossover\_step(self, peptides: List\['Peptide'\], num\_crossovers: int, crossover\_probability\_per\_pair: float) \-\> List\['Peptide'\]**  
    * **Description:** Internal helper to apply crossover to a list of peptides.  
  * **\_execute\_evaluation\_step(self, peptides: List\['Peptide'\], evaluator\_instance: 'Evaluator', round\_identifier: str) \-\> Tuple\[List\['Peptide'\], Optional\[List\[float\]\]\]**  
    * **Description:** Internal helper to execute the evaluation pipeline.  
  * **\_execute\_replenishment\_step(self, peptides: List\['Peptide'\], target\_size: int, gen\_func: Callable) \-\> List\['Peptide'\]**  
    * **Description:** Internal helper to replenish the pool to a target size using the provided generation function.  
  * **\_execute\_truncation\_step(self, peptides: List\['Peptide'\], max\_size: int, scores: Optional\[List\[float\]\]) \-\> List\['Peptide'\]**  
    * **Description:** Internal helper to truncate the pool to a maximum size, optionally using scores for selection.  
  * **run\_round(self, input\_pool\_manager: 'PeptidePoolManager', round\_name: str \= "Evolution Round") \-\> Tuple\['PeptidePoolManager', pd.DataFrame\]**  
    * **Description:** Executes a complete evolutionary round on the peptide pool based on the configured pipeline.  
    * **Arguments:**  
      * input\_pool\_manager (PeptidePoolManager): The PeptidePoolManager of the (n-1) pool.  
      * round\_name (str): Name of this round for logs and evaluator logging.  
    * **Returns:** Tuple\[PeptidePoolManager, pd.DataFrame\]: The resulting (n) pool and a DataFrame of logs for this round.

### **9\. Peptide Generator**

The PeptideGenerator class provides a utility to create random peptide SMILES sequences, useful for initializing peptide pools.

* **class PeptideGenerator**  
  * **\_\_init\_\_(self)**  
    * **Description:** Initializes the generator with the default amino acids from AMINO\_ACID\_DF.  
  * **generate\_random\_peptides(self, num\_peptides: int, min\_length: int, max\_length: int) \-\> List\[List\[str\]\]**  
    * **Description:** Generates a list of random peptide SMILES sequences.  
    * **Arguments:**  
      * num\_peptides (int): The number of peptides to generate.  
      * min\_length (int): The minimum length for generated peptides.  
      * max\_length (int): The maximum length for generated peptides.  
    * **Returns:** List\[List\[str\]\]: A list where each element is a list of SMILES strings representing a peptide sequence.

### **10\. Full Pipeline Examples**

This section provides examples of how to use the various components of PepFlex to construct a complete peptide screening or evolution pipeline.

#### **Example 1: Custom Mutation Function**

This demonstrates how to define and use a custom mutation rule with the PeptideMutator.

* **custom\_reverse\_segment(smiles\_list: List\[str\], start\_idx: int \= 0, end\_idx: Optional\[int\] \= None) \-\> List\[str\]**  
  * **Description:** A custom mutation rule that reverses a segment of the peptide. This function illustrates how a user can define their own mutation logic.  
  * **Arguments:**  
    * smiles\_list (List\[str\]): The SMILES sequence of the peptide.  
    * start\_idx (int): The starting index of the segment to reverse.  
    * end\_idx (Optional\[int\]): The ending index (exclusive) of the segment to reverse. Defaults to the end of the list.  
  * **Returns:** List\[str\]: The peptide SMILES sequence with the specified segment reversed.

#### **Example 2: Evaluator Pipeline Steps**

These functions are examples of steps that can be included in the evaluation\_pipeline of the Evaluator class.

* **add\_length\_feature(df: pd.DataFrame) \-\> pd.DataFrame**  
  * **Description:** Evaluator pipeline step: Adds a 'length' column to the DataFrame based on smiles\_sequence length.  
  * **Arguments:**  
    * df (pd.DataFrame): The input DataFrame.  
  * **Returns:** pd.DataFrame: The DataFrame with the 'length' column added.  
* **add\_dummy\_score\_feature(df: pd.DataFrame) \-\> pd.DataFrame**  
  * **Description:** Evaluator pipeline step: Adds a 'dummy\_score' column for demonstration. This simulates calculating a score based on peptide properties or external models.  
  * **Arguments:**  
    * df (pd.DataFrame): The input DataFrame.  
  * **Returns:** pd.DataFrame: The DataFrame with the 'dummy\_score' column added.  
* **filter\_by\_length\_in\_df(df: pd.DataFrame, min\_len: int \= 7\) \-\> pd.DataFrame**  
  * **Description:** Evaluator pipeline step: Filters the DataFrame to keep peptides above a minimum length.  
  * **Arguments:**  
    * df (pd.DataFrame): The input DataFrame.  
    * min\_len (int): The minimum length to filter by.  
  * **Returns:** pd.DataFrame: The filtered DataFrame.

#### **Example 3: Ranker Function for Evaluator**

This is an example of a ranker\_function that can be passed to the Evaluator class.

* **rank\_by\_dummy\_score\_and\_reconstruct(df: pd.DataFrame, n\_to\_keep: int \= 20\) \-\> Tuple\[List\['Peptide'\], Optional\[List\[float\]\]\]**  
  * **Description:** Ranker function: Ranks peptides by 'dummy\_score' and reconstructs Peptide objects. Returns the top N peptides and their scores.  
  * **Arguments:**  
    * df (pd.DataFrame): The final processed DataFrame from the evaluation pipeline.  
    * n\_to\_keep (int): The number of top peptides to return.  
  * **Returns:** Tuple\[List\['Peptide'\], Optional\[List\[float\]\]\]: A tuple containing:  
    * A list of Peptide objects, ranked by score.  
    * A list of their corresponding scores.

#### **Example 4: Setting up and Running an Evolutionary Pipeline**

This example demonstrates how to combine the PeptideGenerator, PeptideMutator, Evaluator, and PoolRoundProcessor to run a multi-round evolutionary simulation.  
import pandas as pd  
import numpy as np \# Needed for np.random.rand  
import random  
from datetime import datetime  
import copy  
import json  
import uuid  
from typing import List, Tuple, Dict, Callable, Optional, Union

\# Assuming pepflex.py is available and contains all necessary classes and functions.  
\# In a real application, you would import directly from your pepflex module:  
\# from pepflex import (  
\#     AMINO\_ACID\_DF, get\_3L\_from\_smiles, get\_1L\_from\_smiles,  
\#     Peptide, PeptidePoolManager, PeptideMutator, Evaluator, PoolRoundProcessor, PeptideGenerator,  
\#     add\_length\_feature, add\_dummy\_score\_feature, filter\_by\_length\_in\_df,  
\#     rank\_by\_dummy\_score\_and\_reconstruct, peptide\_crossover  
\# )

\# 1\. Initialize Peptide Generator ( You can use your own \- this one if for reference only )  
peptide\_gen \= PeptideGenerator()

\# 2\. Create an initial pool of peptides  
initial\_smiles\_lists \= peptide\_gen.generate\_random\_peptides(num\_peptides=50, min\_length=5, max\_length=15)  
initial\_pool\_manager \= PeptidePoolManager()  
for i, smiles\_list in enumerate(initial\_smiles\_lists):  
    initial\_pool\_manager.add\_peptide(Peptide(smiles\_list, peptide\_id=f"initial\_pep\_{i}"))

print(f"Initial pool size: {initial\_pool\_manager.get\_pool\_size()}")

\# 3\. Configure the Mutator  
mutator \= PeptideMutator()  
mutator.add\_mutation\_rule(mutation\_type='n\_terminal\_addition', probability=0.3)  
mutator.add\_mutation\_rule(mutation\_type='inter\_mutation', probability=0.5)

\# 4\. Configure the Evaluator  
\# Define the evaluation pipeline steps  
evaluation\_pipeline\_steps \= \[  
    add\_length\_feature,  
    add\_dummy\_score\_feature,  
    lambda df: filter\_by\_length\_in\_df(df, min\_len=7) \# Using a lambda for parameter passing  
\]  
\# Define the ranker function  
ranker \= lambda df: rank\_by\_dummy\_score\_and\_reconstruct(df, n\_to\_keep=20)

evaluator \= Evaluator(evaluation\_pipeline=evaluation\_pipeline\_steps, ranker\_function=ranker)

\# 5\. Set up the PoolRoundProcessor  
round\_processor \= PoolRoundProcessor()

\# Set the generation function for replenishment  
round\_processor.set\_generation\_function(  
    lambda num: \[Peptide(s, source\_generation\_params={"type": "replenishment"}) for s in peptide\_gen.generate\_random\_peptides(num, 5, 15)\]  
)

\# Add pipeline steps to the round processor  
\# Note: The step\_function for mutation and crossover needs to be the internal helper,  
\# not the class instance itself, as the helper functions expect \`peptides\` as the first arg.  
\# We pass the mutator instance and other parameters via \*\*step\_params.  
round\_processor.add\_pipeline\_step(  
    step\_type='mutation',  
    step\_function=round\_processor.\_execute\_mutation\_step,  
    name='Apply Mutations',  
    mutator=mutator,  
    probability\_of\_application=0.8 \# Probability that a given peptide will be mutated  
)

round\_processor.add\_pipeline\_step(  
    step\_type='crossover',  
    step\_function=round\_processor.\_execute\_crossover\_step,  
    name='Perform Crossover',  
    num\_crossovers=10, \# Number of crossover events to attempt  
    crossover\_probability\_per\_pair=0.7 \# Probability that a selected pair will actually crossover  
)

round\_processor.add\_pipeline\_step(  
    step\_type='evaluation',  
    step\_function=round\_processor.\_execute\_evaluation\_step,  
    name='Evaluate and Rank Peptides',  
    evaluator\_instance=evaluator  
)

round\_processor.add\_pipeline\_step(  
    step\_type='replenishment',  
    step\_function=round\_processor.\_execute\_replenishment\_step,  
    name='Replenish Pool',  
    target\_size=50 \# Maintain a pool size of 50  
)

round\_processor.add\_pipeline\_step(  
    step\_type='truncation',  
    step\_function=round\_processor.\_execute\_truncation\_step,  
    name='Truncate Pool',  
    max\_size=50 \# Truncate to 50 after all operations  
)

\# 6\. Run multiple rounds of evolution  
current\_pool\_manager \= initial\_pool\_manager  
all\_round\_logs\_df \= pd.DataFrame()

num\_evolution\_rounds \= 3

print("\\n--- Starting Evolutionary Simulation \---")  
for i in range(num\_evolution\_rounds):  
    print(f"\\n===== Running Evolution Round {i+1} \=====")  
    new\_pool\_manager, round\_logs \= round\_processor.run\_round(current\_pool\_manager, round\_name=f"Round\_{i+1}")  
    current\_pool\_manager \= new\_pool\_manager  
    all\_round\_logs\_df \= pd.concat(\[all\_round\_logs\_df, round\_logs\], ignore\_index=True)  
    print(f"End of Round {i+1}. Current pool size: {current\_pool\_manager.get\_pool\_size()}")

print("\\n--- Evolutionary Simulation Finished \---")

\# 7\. Display final results  
print("\\nFinal Top Peptides:")  
final\_peptides \= current\_pool\_manager.get\_all\_peptides()  
for peptide in final\_peptides\[:10\]: \# Display top 10 from final pool  
    print(f"- ID: {peptide.peptide\_id\[:8\]}..., 1L: {peptide.one\_letter\_sequence}, 3L: {peptide.three\_letter\_sequence}, Length: {peptide.length}")

print("\\nEvolutionary Round Logs:")  
print(all\_round\_logs\_df)

**Explanation of the Example:**

1. **Initialization:**  
   * A PeptideGenerator is created to generate random peptide sequences.  
   * An initial\_pool\_manager is populated with 50 randomly generated Peptide objects.  
2. **Mutator Configuration:**  
   * A PeptideMutator is initialized.  
   * Two mutation rules are added: n\_terminal\_addition (adds an amino acid to the N-terminus) and inter\_mutation (substitutes an amino acid within the sequence). Each has a defined probability of occurring.  
3. **Evaluator Configuration:**  
   * The evaluation\_pipeline\_steps list defines the sequence of operations the Evaluator will perform on the peptide DataFrame. Here, it adds a 'length' column, a 'dummy\_score' (simulating a real scoring function), and then filters peptides based on a minimum length.  
   * The ranker function is set to rank\_by\_dummy\_score\_and\_reconstruct, which sorts peptides by their 'dummy\_score' and returns the top 20\.  
   * An Evaluator instance is created with these pipeline steps and the ranker.  
4. **PoolRoundProcessor Setup:**  
   * A PoolRoundProcessor is initialized.  
   * A generation\_function is set, which the replenishment step will use to create new random peptides if the pool size drops below a target.  
   * Various pipeline steps are added to the round\_processor:  
     * **mutation**: Applies the configured mutations from the mutator to the peptides.  
     * **crossover**: Performs crossover operations between peptides to generate new genetic combinations.  
     * **evaluation**: Runs the Evaluator to score and rank peptides. This step is crucial for selection.  
     * **replenishment**: Adds new randomly generated peptides if the pool size falls below target\_size.  
     * **truncation**: Reduces the pool size to max\_size by keeping the highest-scoring peptides.  
5. **Running Evolutionary Rounds:**  
   * The code then enters a loop to run num\_evolution\_rounds (e.g., 3 rounds).  
   * In each round, round\_processor.run\_round() is called with the current peptide pool.  
   * The run\_round method orchestrates the execution of all defined pipeline steps in order.  
   * The current\_pool\_manager is updated with the output of each round, simulating evolution.  
   * Logs from each round are collected.  
6. **Displaying Results:**  
   * Finally, the example prints the details of the top 10 peptides from the final evolved pool and a DataFrame containing the logs from all evolutionary rounds.

This example provides a hands-on illustration of how to set up and execute a basic peptide evolution simulation using PepFlex's modular components. Users can easily modify the mutation rules, evaluation criteria, and round parameters to design more complex and targeted screening pipelines.