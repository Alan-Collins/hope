import sys
import matplotlib.pyplot as plt
from collections import defaultdict


class HomoResult():
	def __init__(self, line):
		bits = line.split()
		self.homo_length = bits[0]
		self.base = bits[1]
		if bits[2] in ["?", "mm", "skip"]:
			self.length = bits[2]
		else:
			self.length = int(bits[2])



def read_hope_file(file):
	data = defaultdict(lambda: defaultdict(list))
	with open(file) as fin:
		for line in fin.readlines()[1:]:
			h = HomoResult(line)
			data[h.homo_length][h.base].append(h.length)

	return data


def main(args=None):
	data = read_hope_file(sys.argv[1])
	to_plot = []
	labels = []
	for length, subd in data.items():
		for base, results in subd.items():
			labels.append(f"{length}_{base}")
			data = []
			for i in results:
				if isinstance(i, int):
					if i > 20:
						data.append(21)
					else:
						data.append(i)
			# to_plot.append([i for i in results if isinstance(i, int)])
			to_plot.append(data)

	# print(labels[0])
	# print(to_plot[0][:10])
	
	# ax.boxplot(to_plot, labels=labels)
	colours = {
	"A": "#CC79A7",
	"T": "#E69F00",
	"C": "#009E73",
	"G": "#0072B2"
	}
	fig, axs = plt.subplots(7,3, figsize=(8, 10), constrained_layout=True)
	i = 0
	j = 0
	for label, data in zip(labels, to_plot):
		ax = axs[i][j]
		if j == 2:
			i += 1
			j = 0
		else:
			j += 1
		homo_length = int(label.split("_")[0])
		homo_base = label.split("_")[1]
		# fig, ax = plt.subplots()
		# fig.set_size_inches(12,8)
		ax.set_title(f"length: {homo_length}, base: {homo_base}")
		ax.hist(
			data,
			bins=[i for i in range(-homo_length-1, 23)],
			log=True,
			align='left',
			color=colours[homo_base]
			)
	plt.savefig(f"{sys.argv[2]}all.png")
	plt.close()

		# homo_length = int(label.split("_")[0])
		# fig, ax = plt.subplots()
		# fig.set_size_inches(12,8)
		# ax.hist(
		# 	data,
		# 	bins=[i for i in range(-homo_length-1,1)],
		# 	log=True,
		# 	align='right'
		# 	)
		# plt.savefig(f"{sys.argv[2]}{label}_deletion.png")
		# plt.close()
	


if __name__ == '__main__':
	main()
