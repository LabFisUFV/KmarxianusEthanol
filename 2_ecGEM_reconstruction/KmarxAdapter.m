classdef KmarxAdapter < ModelAdapter 
	methods
		function obj = KmarxAdapter()
			obj.params.path = fullfile('insert_path_here','KmarxianusEthanol','2_ecGEM_reconstruction');

			obj.params.convGEM = fullfile(obj.params.path,'Kluyveromyces_marxianus-GEM','modelFiles','yml','Kluyveromyces_marxianus-GEM.yml');

			obj.params.sigma = 0.5;

			obj.params.Ptot = 0.546;

			obj.params.f = 0.5;
			
			obj.params.gR_exp = 0.56;

			obj.params.org_name = 'kluyveromyces marxianus';
			
			obj.params.complex.taxonomicID = 4911;

			obj.params.kegg.ID = 'kmx';

			obj.params.kegg.geneID = 'kegg';

			obj.params.uniprot.type = 'proteome';

			obj.params.uniprot.ID = 'UP000065495';

			obj.params.uniprot.geneIDfield = 'gene_orf';

			obj.params.uniprot.reviewed = false;

			obj.params.c_source = 'r_1726'; 

			obj.params.bioRxn = 'r_1912';

			obj.params.enzyme_comp = 'Cytoplasm';			
		end
	
		function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)
			spont = contains(model.rxnNames,'spontaneous');
			spontRxnNames = model.rxnNames(spont);
		end
		
	end
end
