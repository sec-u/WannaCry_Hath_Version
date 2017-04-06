
import networkx as nx
import re
import copy

class ExtractSynthesisComponents(object):
	def _find_components(self,r,met_type,met_type2):
		rxns=[]
		mets=[]
		if (set(self.rxns[r]['mets'][met_type])<set(self.internalcompounds)) is True:
			rxns.append(r)
			mets= self.rxns[r]['mets'][met_type2]
		return(mets,rxns)
	
	def extract_synth_rxn_met(self):
		countiter=0
		loop=1
		e_rxns=self.externalrxns
		e_rxnscopy=copy.deepcopy(e_rxns)

		while loop == 1:
			countiter+=1
			rxns_temp=[]
			mets_temp=[]
			for r in e_rxns:
				if self.rxns[r]['reversibility']=='false':
					list_mets,rxn=self._find_components(r,'reactants','products')
					mets_temp=mets_temp+list_mets
					rxns_temp=rxns_temp+rxn
				else:
					"""FORWARD DIRECTION"""
					list_mets,rxn=self._find_components(r,'reactants','products')			
					mets_temp=mets_temp+list_mets
					rxns_temp=rxns_temp+rxn
					"""REVERSE DIRECTION"""
					list_mets,rxn=self._find_components(r,'products','reactants')
					mets_temp=mets_temp+list_mets
					rxns_temp=rxns_temp+rxn

			self.internalcompounds=self.internalcompounds+mets_temp
			self.internalcompounds=set(self.internalcompounds) ##REMOVES DUPLICATES
			self.internalcompounds=list(self.internalcompounds) ##CONVERTS BACK TO LIST
			
			rxns_temp=set(rxns_temp) ##REMOVES DUPLICATES( MAY HAVE BEEN ADDED IF RXN IS REVERSIBLE)
			rxns_temp=list(rxns_temp) ##CONVERTS BACK TO LIST
			
			for e in rxns_temp:
				e_rxns.remove(e)

			self.rxns_iter[countiter]=rxns_temp
			self.mets_prod_iter[countiter]=mets_temp
			if self.out_met in mets_temp:
				print 'Iteration stopped because target metabolite could be produced'
				loop=0
			elif len(self.mets_prod_iter) >1 and len(self.rxns_iter)>1:
				if len(self.rxns_iter[countiter]) >= len(self.rxns_iter[countiter-1]):
					print str(len(self.rxns_iter[countiter]))+'\t'+str(len(self.rxns_iter[countiter-1]))
					print 'Iteration stopped because rxn list not getting smaller, target metabolite could not be identified'
					loop=0

		print 'Exited loop on iteration '+str(countiter)
		print self.rxns_iter
		print self.mets_prod_iter

	def __init__(self,er,out_met):
		self.out_met=out_met
		self.externalrxns=er.externalrxns
		self.internalcompounds=er.internalcompounds
		self.rxns=er.rxns
		self.rxns_iter={}
		self.mets_prod_iter={}
		self.extract_synth_rxn_met()
