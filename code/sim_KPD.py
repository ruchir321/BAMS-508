import numpy as np
import pandas as pd
from scipy.optimize import milp, LinearConstraint, Bounds
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Tuple, List, Dict
import time

class KindneyDonationSimulator:
    def __init__(self):
         # US population blood type distribution (including Rh factor)
        self.blood_type_probs = {
            'AB': 0.034, 'A': 0.357, 'B': 0.085, 'O': 0.374, 
            'AB-': 0.006, 'A-': 0.063, 'B-': 0.015, 'O-': 0.066
        }

        self.blood_types = list(self.blood_type_probs.keys())
        self.blood_probs = list(self.blood_type_probs.values())
        
        # blood compatibility matrix in order: [AB, A, B, O, AB-, A-, B- , O-]
        self.M_blood = np.array(
               [[1,0,0,0,0,0,0,0],
                [1,1,0,0,0,0,0,0],
                [1,0,1,0,0,0,0,0],
                [1,1,1,1,0,0,0,0],
                [1,0,0,0,1,0,0,0],
                [1,0,0,0,1,1,0,0],
                [1,0,0,0,1,0,1,0],
                [1,1,1,1,1,1,1,1]], dtype=int)
        
        self.M_tissue = np.identity(n=5,dtype=int) # tissue compatibility matrix simplified to 5 tissue types

        self.tissue_probs = [0.2, 0.2, 0.2, 0.2, 0.2]

    def generate_random_data(self, D: int, P: int, seed: int=None) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Generates data from realistic distributions for columns: DonorBloodType, DonorTissueType, PatientBloodType, PatientTissueType
        """
        
        if seed is not None:
            np.random.seed(seed)
        
        donor_blood_indices = np.random.choice(
            a=len(self.blood_types),
            size=D,
            p=self.blood_probs
        )

        patient_blood_indices = np.random.choice(
            a=len(self.blood_types),
            size=P,
            p=self.blood_probs
        )

        donor_tissue_indices = np.random.choice(
            a=5,
            size=D,
            p=self.tissue_probs
        )

        patient_tissue_indices = np.random.choice(
            a=5,
            size=P,
            p=self.tissue_probs
        )

        # One Hot encoding
        bt_D = np.zeros((D, 8), dtype=int)
        bt_P = np.zeros((P, 8), dtype=int)

        t_D = np.zeros((D, 5), dtype=int)
        t_P = np.zeros((P, 5), dtype=int)

        for i, blood_idx in enumerate(donor_blood_indices):
            bt_D[i, blood_idx] = 1
        
        for i, blood_idx in enumerate(patient_blood_indices):
            bt_P[i, blood_idx] = 1

        for i, tissue_idx in enumerate(donor_tissue_indices):
            t_D[i, tissue_idx] = 1

        for i, tissue_idx in enumerate(patient_tissue_indices):
            t_P[i, tissue_idx] = 1

        return bt_D, bt_P, t_D, t_P
    
    def optimize(self, bt_D: np.ndarray, bt_P: np.ndarray, t_D: np.ndarray, t_P: np.ndarray) -> Dict:
        """
        Use Integer linear programming to find maximum compatible transplants possible
        """

        pairs = (bt_D @ self.M_blood @ bt_P.T) * (t_D @ self.M_tissue @ t_P.T)

        D = pairs.shape[0]
        P = pairs.shape[1]

        coefficients = -pairs.flatten()

        A_donor = np.zeros(shape=(D, D*P), dtype=int)
        A_patient = np.zeros(shape=(P, D*P), dtype=int)

        for i in range(D):
            for j in range(P):
                A_donor[i, i * P + j] = 1

        for j in range(P):
            for i in range(D):
                A_patient[j, i * P + j] = 1

        A = np.vstack((A_donor, A_patient))
        b_ub = np.ones(D + P)

        constraints = LinearConstraint(
            A=A,
            lb=0,
            ub=b_ub)

        bounds = Bounds(
            lb=0,
            ub=1)
        
        integrality = np.ones(D*P, dtype=int)
        
        result_milp = milp(
            c=coefficients,
            bounds=bounds,
            constraints=constraints,
            integrality=integrality)
        

        selection_matrix = result_milp.x.reshape(D, P)
        n_matches = int(np.sum(selection_matrix))

        matches = []
        for i in range(D):
            for j in range(P):
                if selection_matrix[i, j] == 1:
                    matches.append({
                        'donor_idx': i,
                        'patient_idx': j,
                        'compatibility_score': pairs[i, j]
                    })
        
        return {
            "success": result_milp.success,
            "n_matches": n_matches,
            "selection_matrix": selection_matrix,
            "pairs_matrix": pairs,
            "total_score": -result_milp.fun if result_milp.success else 0
        }
    

    def run_simulation(self, donor_counts: List[int], patient_counts: List[int], n_trials: int=10) -> pd.DataFrame:
        """
        Run simulation across different donor/patient ratios
        """

        results = []

        for n_donors in donor_counts:
            for n_patients in patient_counts:
                for trial in range(n_trials):
                    seed = 42*trial + 100*n_donors + n_patients
                    bt_D, bt_P, t_D, t_P = self.generate_random_data(D=n_donors, P=n_patients, seed=seed)

                    result = self.optimize(bt_D, bt_P, t_D, t_P)

                    if result["success"]:
                        results.append({
                            "n_donors": n_donors,
                            "n_patients": n_patients,
                            "trial_number": trial,
                            "n_matches": result["n_matches"],
                            "total_score": result["total_score"],
                            "match_rate_donors": result["n_matches"] / n_donors,
                            "match_rate_patients": result["n_matches"] / n_patients,
                            "donor_patient_ratio": n_donors / n_patients
                            })
        
        return pd.DataFrame(results)
    


sim = KindneyDonationSimulator()

print("Running simulation...")

donor_counts = (np.linspace(start=50,stop=100,num=11, dtype=int)).tolist()
patient_counts = np.linspace(start=50,stop=100, num=11, dtype=int).tolist()

# donor_counts = [10, 150, 200, 25, 30, 35, 40]
# patient_counts = [10, 150, 20, 250, 30, 35, 40]


start_time = time.time()
sim_result = sim.run_simulation(donor_counts=donor_counts, patient_counts=patient_counts, n_trials=5)
end_time = time.time()
run_time = end_time - start_time


print(f"Simulation complete: {len(sim_result)} experiments in {run_time:0.2f} seconds")


plt.figure(figsize=(18,12))

# 1.How does changing supply/demand ratio affect matches
plt.subplot(2,3,1)
avg_results = sim_result.groupby(by=["n_donors", "n_patients"]).agg(
    {
        "match_rate_donors": "mean",
        "match_rate_patients": "mean",
        "donor_patient_ratio": "mean",
        "n_matches": "mean",
    }
).reset_index()

plt.scatter(avg_results["donor_patient_ratio"], avg_results["match_rate_donors"], label="Donor Match Ratio", alpha=0.8)
plt.scatter(avg_results["donor_patient_ratio"], avg_results["match_rate_patients"], label="Patient Match Ratio", alpha=0.8)

plt.xlabel("Donor Patient Ratio")
plt.ylabel("Match Rate")

plt.title("Match Rate VS Donor Patient Ratio")
plt.legend()
plt.grid(True, alpha=0.3)

# 2. How does changing the problem size affect matches
plt.subplot(2, 3, 2)

avg_results['total_people'] = avg_results['n_donors'] + avg_results['n_patients']

plt.scatter(avg_results['total_people'], avg_results['n_matches'], c=avg_results['donor_patient_ratio'], cmap='viridis', alpha=0.7)

plt.colorbar(label='Donor/Patient Ratio')
plt.xlabel('Total People (Donors + Patients)')
plt.ylabel('Number of Matches')
plt.title('Matches vs Problem Size')
plt.grid(True, alpha=0.3)

# Plot 3: Heatmap of match rates
plt.subplot(2, 3, 3)

pivot_matches = avg_results.pivot(index='n_patients', columns='n_donors', values='n_matches')

sns.heatmap(pivot_matches, annot=True, fmt='.1f', cmap='YlOrRd')
plt.title('Average Number of Matches')
plt.ylabel('Number of Patients')
plt.xlabel('Number of Donors')

# Create text summary
plt.subplot(2, 3, 5)
plt.text(0.1, 0.8, f"Simulation Summary:\n" +
         f"Total experiments: {len(sim_result)}\n" +
         f"Avg matches: {sim_result['n_matches'].mean():.1f} ±{sim_result['n_matches'].std():.1f}\n" +
         f"Avg donor match rate: {sim_result['match_rate_donors'].mean():.3f}\n" +
         f"Avg patient match rate: {sim_result['match_rate_patients'].mean():.3f}\n\n" +
         f"Best scenario:\n" +
         f"Max matches: {sim_result['n_matches'].max()}\n" +
         f"At ratio: {sim_result.loc[sim_result['n_matches'].idxmax(), 'donor_patient_ratio']:.2f}",
         transform=plt.gca().transAxes, fontsize=10, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
plt.axis('off')

plt.tight_layout()
plt.savefig('output/simulation_results.png')
plt.show()

print("\n2. Optimal Scenarios (Top 10 by total matches):")
top_scenarios = sim_result.nlargest(10, 'n_matches')[
    ['n_donors', 'n_patients', 'n_matches', 'donor_patient_ratio', 'match_rate_donors', 'match_rate_patients']
]
print(top_scenarios)