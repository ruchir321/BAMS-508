# Kidney Paired Donation (KPD) Optimization Analysis

This repository contains a simulation and optimization model for the Kidney Paired Donation (KPD) program. The project analyzes the dynamics of kidney exchanges to identify scenarios that maximize the number of successful transplants. The underlying model, detailed in `code/sim_KPD.py`, uses integer linear programming to resolve compatible matches based on blood and tissue type compatibility.

## Background

The Kidney Paired Donation (KPD) program, managed by Canadian Blood Services (CBS), facilitates a voluntary kidney exchange for incompatible donor-patient pairs. When a patient cannot receive a kidney from their intended donor due to biological incompatibility, KPD provides an opportunity to find a compatible match with another pair in the system. This project models a simplified version of this process, focusing on the core optimization problem of maximizing transplants through exchange cycles. Compatibility is determined by both blood type and a direct match in tissue type.

## How to Run

1.  **Set up the environment:**
    ```bash
    pip install -r requirements.txt
    ```
2.  **Run the simulation:**
    ```bash
    python code/sim_KDP.py
    ```

## Output and Analysis

The simulation generates a comprehensive analysis of matching outcomes under varying conditions, summarized in `output/simulation_results.png`.

The key findings from the 605 simulation experiments are:
- The total number of transplants is strongly correlated with the size of the participant pool (both donors and patients).
- The ratio of donors to patients significantly impacts the efficiency of the matching process.

**Analysis of Results:**

The simulation results, visualized in the output plots, provide critical insights for optimizing the KDP pool.

-   **Match Rate Dynamics:** The 'Match Rate VS Donor Patient Ratio' plot illustrates a fundamental trade-off. As the donor-to-patient ratio increases, the likelihood of a patient finding a match (Patient Match Rate) improves, while the likelihood of a given donor being used in a transplant (Donor Match Rate) decreases.
-   **Optimal Pool Composition:** The 'Average Number of Matches' heatmap and the 'Simulation Summary' provide a clear recommendation. The maximum number of transplants is achieved when the number of donors is close to the number of patients. The simulation identified an optimal **donor-to-patient ratio of 0.95**, which yielded a maximum of 87 successful matches.

**Recommendation:**

To maximize the number of successful transplants, the KDP program should strategically manage its pool of participants to maintain a donor-to-patient ratio as close to 1.0 as possible, with a slight surplus of patients being beneficial. The analysis indicates that a ratio of **0.95 donors to patients** represents the most efficient scenario for maximizing successful exchanges. This ensures a large and balanced pool, increasing the probability of finding compatible matches for the greatest number of patients.
