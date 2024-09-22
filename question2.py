import matplotlib.pyplot as plt

# Hardcoded time values for each algorithm (in seconds, using double)
Brute_Force = 12.5   # Example time for Brute Force in seconds
Median_String = 8.9  # Example time for Median String in seconds
Gibbs_Sampling = 15.7  # Example time for Gibbs Sampling (> 500 iterations) in seconds

# Method names and corresponding times
methods = ['Brute Force', 'Median String', 'Gibbs Sampling']
times = [Brute_Force, Median_String, Gibbs_Sampling]

# Plotting the time comparison
plt.figure(figsize=(8, 6))
plt.bar(methods, times, color=['blue', 'green', 'red'])
plt.title('Time Comparison of Motif Finding Algorithms')
plt.xlabel('Algorithm')
plt.ylabel('Time (seconds)')
plt.show()
