import csv

# read allGoodValues.dat to a list of lists
datContent = [i.strip().split() for i in open("./allGoodValues.dat").readlines()]

# write it as a new CSV file
with open("./allGoodValues.csv", "wb") as f:
    writer = csv.writer(f)
    writer.writerows(datContent)