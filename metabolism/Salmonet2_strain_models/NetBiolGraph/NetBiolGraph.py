#!/usr/bin/python
import json
import _pickle as cPickle

__all__ = ["NetBiolGraph", ]

class NetBiolGraph(dict):

	def __init__(self):
		self._dict = {}
		self._edge_num = 0

	def __str__(self):
		return "NetBiolGraph object, Nodes: %s, Edges: %s" % (len(self), self._edge_num)

	def n_node(self, name, **kwargs):
		if name in self:
			return self[name]
		self[name] = {}
		node = self[name]
		node["name"] = name
		node["aliases"] = {}
		node["species"] = 0
		node["label"] = name
		node["type"] = "protein"
		node["is_tf"] = False
		node["properties"] = {}
		node["edges"] = []
		for key, value in kwargs.items():
			node[key] = value
		return node

	def n_edge(self, source="", target="", **kwargs):
		edge = {}
		edge["name"] = target
		edge["label"] = "%s - %s" % (source, target)
		edge["direction"] = 0
		edge["xrefs"] = {"db": [], "pubmed": []}
		edge["directness"] = 0
		edge["effect"] = 0
		edge["scores"] = []
		for key, value in kwargs.items():
			edge[key] = value
		return edge

	def add_edge(self, source, target, **kwargs):
		source_node = self.n_node(source)
		target_node = self.n_node(target)
		edge_a = self.n_edge(source, target, **kwargs)
		edge_b = self.n_edge(target, source, **kwargs)
		source_node["edges"].append(edge_a)
		target_node["edges"].append(edge_b)
		self._edge_num += 1
		return None

	def edit_edge(self, source, target, **kwargs):
		source_node = self.n_node(source)
		target_node = self.n_node(target)
		edge_a = None
		edge_b = None
		for e in source_node["edges"]:
			if e["name"] == target:
				edge_a = e
		for e in target_node["edges"]:
			if e["name"] == source:
				edge_b = e
		if edge_a is None or edge_b is None:
			raise ValueError()
		for key, value in kwargs.items():
			if key == "direction" and value != 0:
				edge_a[key] = value
				edge_b[key] = 0-value
			elif type(value) == type([]):
				if key not in edge_a:
					edge_a[key] = [] 
				if key not in edge_b:
					edge_b[key] = []
				for v in value:
					edge_a[key].append(v)
					edge_b[key].append(v)
			elif type(value) == type({}):
				if key not in edge_a:
					edge_a[key] = {}
				if key not in edge_b:
					edge_b[key] = {}
				for k, v in value.items():
					edge_a[key][k] = v
					edge_b[key][k] = v
			else:
				edge_a[key] = value
				edge_b[key] = value

	def write_file(self, file_name, file_type="json", json_indent=4, pickle_protocol=0):
		if file_type == "json":
			with open(file_name, "w") as f:
				f.write(json.dumps(self, indent=json_indent))
		elif file_type == "pickle":
			with open(file_name, "wb") as f:
				pickle.dump(self, f, pickle_protocol)
		return True

	@classmethod
	def read_file(klass, file_name, file_type="json"):
		if file_type == "json":
			return json.load(open(file_name))
		elif file_type == "pickle":
			return pickle.load(open(file_name, "rb"))