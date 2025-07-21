
# PepFlex: A Modular Framework for Peptide Screening and Evolution

PepFlex is a powerful Python module designed to simplify *in silico* peptide screening and evolutionary studies. It provides a flexible and extensible framework that empowers researchers to manage amino acid data, assemble intricate peptide structures, apply diverse genetic operations, and guide evolutionary processes through a highly customizable evaluation pipeline. Built with modularity in mind, PepFlex enables the construction of sophisticated workflows, allowing for rapid prototyping and detailed exploration of peptide landscapes.

## Guiding Evolution: The Round Processor

At the heart of PepFlex's evolutionary capabilities lies **the PoolRoundProcessor**. This sophisticated component acts as the orchestrator for a single generation of peptide evolution, seamlessly guiding a population of peptides through a series of transformative steps. It's like the conductor of an orchestra, ensuring each genetic operation and selection pressure plays its part in shaping the next generation. The true power of the PoolRoundProcessor lies in its completely configurable nature; you, the user, define the exact sequence and even the repetition of events within each round. This allows for tailored evolutionary pathways that mimic natural selection or explore specific design spaces with unparalleled precision.

When you initiate a round of evolution, the PoolRoundProcessor takes the current population of peptides, encapsulated within a **PeptidePoolManager instance**, as its input. After meticulously applying all the defined steps : mutation, crossover, evaluation, replenishment, and truncation; it then yields a new, evolved PeptidePoolManager instance as its primary output, representing the next generation of peptides. **Additionally, it provides a detailed log of the round's events**, offering valuable insights into the evolutionary progress.

<img width="400" height="400" alt="poolroundprocessor" src="https://github.com/user-attachments/assets/2e2815b9-1120-4335-a662-32c604cef9e1" />


### Shaping Diversity: Mutation

The **Mutation** step introduces genetic variation into the peptide population, mirroring the random changes that occur in nature. Here, individual amino acids within a peptide sequence can be altered, inserted, or deleted. You can define specific types of mutations, such as changing one amino acid to another, adding new amino acids to either end of a peptide, or even more complex intra-sequence rearrangements. This crucial step ensures that the peptide pool constantly explores new chemical possibilities, preventing stagnation and fostering innovation in the search for desired properties. You can strategically place mutation steps at any point in your round, or even apply multiple, distinct mutation stages.

### Blending Traits: Crossover

The **Crossover** step allows for the recombination of genetic material between different parent peptides. Imagine two promising peptides, each with valuable characteristics; crossover enables segments of their sequences to be exchanged, creating entirely new "child" peptides that inherit traits from both parents. This process is vital for exploring diverse combinations of existing features, accelerating the discovery of novel and potent peptide sequences by blending successful elements from different lineages. Like all steps, you decide where crossover fits into your evolutionary strategy within a round.

### Measuring Success: Evaluation

The **Evaluation** phase is where the "fitness" of each peptide is assessed. This is the selective pressure that drives the evolutionary process towards your desired outcomes. At its core, the **Evaluator is a powerful, customizable pipeline** composed of a series of user-defined functions (often referred to as E0 to En steps). The process begins by converting the population of peptides into a structured dataset. Then, each function in your defined pipeline takes this dataset as input and can perform virtually any operation: from calculating simple features like length or molecular weight, to running complex machine learning models predicting activity or binding affinity, to filtering out peptides that don't meet certain criteria. Each step in this pipeline transforms or enriches the data, leading to a comprehensive understanding of each peptide's potential. Finally, a dedicated "ranker" function uses this refined information to score and order the peptides, identifying the most promising candidates. This modularity means you have complete freedom to define precisely how "success" is measured and what characteristics are prioritized in your evolutionary search. You might even incorporate multiple evaluation steps throughout a round, perhaps a quick preliminary screen followed by a more in-depth analysis for promising candidates.

<img width="594" height="177" alt="evaluator structure" src="https://github.com/user-attachments/assets/ea55cf6a-30b8-45ad-8139-cd751fc00fcd" />

### Sustaining the Pool: Replenishment

To maintain a healthy and vibrant evolutionary population, the **Replenishment** step ensures that the peptide pool never dwindles too low. If the population size falls below a desired threshold, this mechanism introduces new peptides into the mix. Crucially, **the way these new peptides are generated is entirely up to you.** While a basic random generator is provided, you can define your own custom function to create new peptides. This flexibility allows for targeted replenishment strategies, perhaps introducing peptides with specific motifs, or drawing from external databases, ensuring the continued exploration of diverse and relevant peptide chemistries. This influx of fresh sequences prevents premature convergence, broadens the genetic diversity, and ensures that the evolutionary process continues to explore novel regions of the peptide space.

### Refining the Elite: Truncation

The **Truncation** step is where the principle of "survival of the fittest" comes into play. After other operations have taken place, the population is sorted based on their assessed "fitness scores." Truncation then prunes the pool, removing the least promising peptides and retaining only the top performers. This ensures that the limited resources of the simulation are focused on the most promising candidates, concentrating the evolutionary pressure and effectively driving the population towards peptides with superior properties. Like evaluation, you can use truncation multiple times within a round, perhaps to aggressively filter early on, then to perform a final selection at the end.

### What PepFlex Brings to the Researcher: Empowering Discovery, Not Rebuilding

PepFlex is designed to be an enabling platform for discovery, not a toolkit that forces you to rebuild from scratch. For researchers, it dramatically lowers the barrier to entry for conducting sophisticated in silico peptide evolutionary studies. Instead of engineering complex simulation infrastructures, managing data flow between disparate scripts, and developing every genetic operator or evaluation metric from the ground up, PepFlex provides a robust, pre-built backbone.

**So its modular architecture means you get the best of both worlds**:

First, Ready-to-Use Components: PepFlex comes with essential components like a peptide generator, common mutation types, a pool manager, and a powerful evaluator framework already established. This means you can get a basic evolutionary simulation up and running with minimal effort and without writing hundreds of lines of boilerplate code.

Next, Points of Customization: When your research demands unique approaches, PepFlex offers specific, well-defined points where you can inject your own custom logic. Whether it's a novel scoring function, an innovative mutation rule, or a specialized way to generate new peptides, you only need to build that specific piece, not the entire evolutionary engine.

This approach means you can rapidly prototype and test different hypotheses about peptide design, explore vast chemical spaces efficiently, and tailor your search to very specific therapeutic or industrial targets. PepFlex empowers you to focus directly on your scientific questions, dedicating your valuable time to designing insightful experiments and interpreting results, rather than to the tedious work of framework development. It provides the solid foundation you need, coupled with the freedom to innovate where it truly matters.

### Example usage  : Construction of a pipeline

First, you can install pepflex with the pip package manager : 
```
pip install pepflex
```
The following example ( **From DOCUMENTATION.md** ) demonstrates how to combine the PeptideGenerator, PeptideMutator, Evaluator, and PoolRoundProcessor to run a multi-round evolutionary simulation. .  

```
import pandas as pd
import numpy as np  # Needed for np.random.rand
import random
from datetime import datetime
import copy
import json
import uuid
from typing import List, Tuple, Dict, Callable, Optional, Union

# In a real application, you would import directly from your pepflex module:
# from pepflex import (
#     AMINO_ACID_DF, get_3L_from_smiles, get_1L_from_smiles,
#     Peptide, PeptidePoolManager, PeptideMutator, Evaluator, PoolRoundProcessor, PeptideGenerator,
#     add_length_feature, add_dummy_score_feature, filter_by_length_in_df,
#     rank_by_dummy_score_and_reconstruct, peptide_crossover
# )

# 1. Initialize Peptide Generator ( You can use your own - this one is for reference only )
peptide_gen = PeptideGenerator()

# 2. Create an initial pool of peptides
initial_smiles_lists = peptide_gen.generate_random_peptides(num_peptides=50, min_length=5, max_length=15)
initial_pool_manager = PeptidePoolManager()
for i, smiles_list in enumerate(initial_smiles_lists):
    initial_pool_manager.add_peptide(Peptide(smiles_list, peptide_id=f"initial_pep_{i}"))

print(f"Initial pool size: {initial_pool_manager.get_pool_size()}")

# 3. Configure the Mutator
mutator = PeptideMutator()
mutator.add_mutation_rule(mutation_type='n_terminal_addition', probability=0.3)
mutator.add_mutation_rule(mutation_type='inter_mutation', probability=0.5)

# 4. Configure the Evaluator
# Define the evaluation pipeline steps
evaluation_pipeline_steps = [
    add_length_feature,
    add_dummy_score_feature,
    lambda df: filter_by_length_in_df(df, min_len=7)  # Using a lambda for parameter passing
]
# Define the ranker function
ranker = lambda df: rank_by_dummy_score_and_reconstruct(df, n_to_keep=20)

evaluator = Evaluator(evaluation_pipeline=evaluation_pipeline_steps, ranker_function=ranker)

# 5. Set up the PoolRoundProcessor
round_processor = PoolRoundProcessor()

# Set the generation function for replenishment
round_processor.set_generation_function(
    lambda num: [Peptide(s, source_generation_params={"type": "replenishment"})
                 for s in peptide_gen.generate_random_peptides(num, 5, 15)]
)

# Add pipeline steps to the round processor
round_processor.add_pipeline_step(
    step_type='mutation',
    step_function=round_processor._execute_mutation_step,
    name='Apply Mutations',
    mutator=mutator,
    probability_of_application=0.8  # Probability that a given peptide will be mutated
)

round_processor.add_pipeline_step(
    step_type='crossover',
    step_function=round_processor._execute_crossover_step,
    name='Perform Crossover',
    num_crossovers=10,  # Number of crossover events to attempt
    crossover_probability_per_pair=0.7  # Probability that a selected pair will actually crossover
)

round_processor.add_pipeline_step(
    step_type='evaluation',
    step_function=round_processor._execute_evaluation_step,
    name='Evaluate and Rank Peptides',
    evaluator_instance=evaluator
)

round_processor.add_pipeline_step(
    step_type='replenishment',
    step_function=round_processor._execute_replenishment_step,
    name='Replenish Pool',
    target_size=50  # Maintain a pool size of 50
)

round_processor.add_pipeline_step(
    step_type='truncation',
    step_function=round_processor._execute_truncation_step,
    name='Truncate Pool',
    max_size=50  # Truncate to 50 after all operations
)

# 6. Run multiple rounds of evolution
current_pool_manager = initial_pool_manager
all_round_logs_df = pd.DataFrame()

num_evolution_rounds = 3

print("\n--- Starting Evolutionary Simulation ---")
for i in range(num_evolution_rounds):
    print(f"\n===== Running Evolution Round {i+1} =====")
    new_pool_manager, round_logs = round_processor.run_round(current_pool_manager, round_name=f"Round_{i+1}")
    current_pool_manager = new_pool_manager
    all_round_logs_df = pd.concat([all_round_logs_df, round_logs], ignore_index=True)
    print(f"End of Round {i+1}. Current pool size: {current_pool_manager.get_pool_size()}")

print("\n--- Evolutionary Simulation Finished ---")

# 7. Display final results
print("\nFinal Top Peptides:")
final_peptides = current_pool_manager.get_all_peptides()
for peptide in final_peptides[:10]:  # Display top 10 from final pool
    print(f"- ID: {peptide.peptide_id[:8]}..., 1L: {peptide.one_letter_sequence}, "
          f"3L: {peptide.three_letter_sequence}, Length: {peptide.length}")

print("\nEvolutionary Round Logs:")
print(all_round_logs_df)

```
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

