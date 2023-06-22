import yaml


class Config(object):
	def __init__(self, fname):
		self.fname = fname
		self.read_data()

	def read_data(self):

		with open(self.fname, "r") as f:
		    self.v = yaml.load(f)
