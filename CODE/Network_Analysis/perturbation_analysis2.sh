whoami
date
start=`date +%s`


# declare conditions to analyze
declare -a cell_type=("case_vs_controls")


for cell in "${cell_type[@]}"; do

echo -e "########################################################\nPerturbation analysis for $cell\n########################################################"

# prepare input temporary names
temp_nodemap="${cell}_network_nodes.txt"
temp_interact="${cell}_network_interactions.txt"
temp_gene_list="${cell}_genes_list.txt"
temp_gene_expression="${cell}_genes_expression.txt"
# prepare output temporary names
temp_adjacency="${cell}_adjacency.txt"
temp_networkphenotype1="${cell}_NetworkPhenotype1.txt"
temp_networkphenotype2="${cell}_NetworkPhenotype2.txt"
temp_CommonNetworkGenerator="${cell}_CommonNetworkGenerator_Output.txt"
temp_DiffNetworkGenerator="${cell}_DifferentialNetworkGenerator_Output.txt"
temp_pos="${cell}_pos.txt"
temp_neg="${cell}_neg.txt"
temp_SteadyStateCalculatorN1="${cell}_SteadyStateCalculatorN1.txt"
temp_SteadyStateCalculatorN2="${cell}_SteadyStateCalculatorN2.txt"
temp_PerturbagenListGeneratorN1="${cell}_PerturbagenListGeneratorN1.txt"
temp_PerturbagenListGeneratorN2="${cell}_PerturbagenListGeneratorN2.txt"
temp_BruteForcePerturbationsUpdatedN1_1="${cell}_BruteForcePerturbationsUpdatedN1_1.txt"
temp_BruteForcePerturbationsUpdatedN1_2="${cell}_BruteForcePerturbationsUpdatedN1_2.txt"
temp_BruteForcePerturbationsUpdatedN1_3="${cell}_BruteForcePerturbationsUpdatedN1_3.txt"

echo $"$temp_nodemap"

## Pre-process input data
java -jar Preprocessor.jar . "$temp_nodemap" "$temp_interact" "$temp_gene_list"
mv "adjacency.txt" "$temp_adjacency"

## run differential network analysis
java -jar DifferentialNetworkAnalysis.jar "$temp_gene_expression" "$temp_adjacency"  GAResult.txt 0 true 1000 50 .

## rename NetworkPhenotype1.txt & NetworkPhenotype2.txt to case_vs_controls_NetworkPhenotype1.txt & case_vs_controls_NetworkPhenotype2.txt
mv NetworkPhenotype1.txt "$temp_networkphenotype1"
mv NetworkPhenotype2.txt "$temp_networkphenotype2"

### Create here Cytoscape vizualisation if needed

## Perturbation analysis
java -jar CommonNetworkGenerator.jar "$temp_networkphenotype1" "$temp_networkphenotype2" "$temp_CommonNetworkGenerator"
java -jar DifferentialNetworkGenerator.jar "$temp_networkphenotype1" "$temp_networkphenotype2" "$temp_DiffNetworkGenerator"
java -jar ComputeCycles.jar "$temp_CommonNetworkGenerator" "$temp_gene_expression" "$temp_pos" "$temp_neg"

## check that there are positive or negative cycles (if not, the perturbation analysis may not be feasible, e.g. if the network is too small --> consider to relax p-value cut-off)

# cat case_vs_controls_pos.txt
# cat case_vs_controls_neg.txt

java -jar SteadyStateCalculator.jar "$temp_gene_expression" "$temp_networkphenotype1" 1 "$temp_SteadyStateCalculatorN1"
java -jar SteadyStateCalculator.jar "$temp_gene_expression" "$temp_networkphenotype1" 2 "$temp_SteadyStateCalculatorN2"

java -jar PerturbagenListGenerator.jar "$temp_pos" "$temp_neg" "$temp_DiffNetworkGenerator" "$temp_SteadyStateCalculatorN1" "$temp_PerturbagenListGeneratorN1"
java -jar PerturbagenListGenerator.jar "$temp_pos" "$temp_neg" "$temp_DiffNetworkGenerator" "$temp_SteadyStateCalculatorN2" "$temp_PerturbagenListGeneratorN2"


java -jar BruteForcePerturbationsUpdated.jar "$temp_gene_expression" "$temp_networkphenotype1" 1 "$temp_PerturbagenListGeneratorN1" 1 500000 "$temp_BruteForcePerturbationsUpdatedN1_1"
java -jar BruteForcePerturbationsUpdated.jar "$temp_gene_expression" "$temp_networkphenotype1" 1 "$temp_PerturbagenListGeneratorN1" 2 500000 "$temp_BruteForcePerturbationsUpdatedN1_2"
java -jar BruteForcePerturbationsUpdated.jar "$temp_gene_expression" "$temp_networkphenotype1" 1 "$temp_PerturbagenListGeneratorN1" 3 500000 "$temp_BruteForcePerturbationsUpdatedN1_3"

# sort -rn "$temp_BruteForcePerturbationsUpdatedN1_1" > Z.txt
# mv Z.txt "$temp_BruteForcePerturbationsUpdatedN1_1"


sort -rn "$temp_BruteForcePerturbationsUpdatedN1_1" > Z.txt
mv Z.txt "$temp_BruteForcePerturbationsUpdatedN1_1"

sort -rn "$temp_BruteForcePerturbationsUpdatedN1_2" > Z.txt
mv Z.txt "$temp_BruteForcePerturbationsUpdatedN1_2"

sort -rn "$temp_BruteForcePerturbationsUpdatedN1_3" > Z.txt
mv Z.txt "$temp_BruteForcePerturbationsUpdatedN1_3"


done

echo "Analysis is done"
# date

end=`date +%s`
echo ########################################################
echo Execution time was `expr $end - $start` seconds.
echo ########################################################
