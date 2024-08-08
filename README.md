# Early Onset IMC
Colorectal cancer (CRC) incidence in patients under 50 has been increasing steadily for the past 20 years. From 2000 to 2013, CRC incidence in young patients (< 50 years old) has increased more than 20%. Most distrurbingly, early-onset CRC (EOCRC) cases are frequently detected at an advanced stage, stage II and III cancer account for about 70% of all EOCRCs. As a consequence, younger patients tend to have poorer outcomes. 
Here, imaging mass cytometry (IMC) was emplyed. IMC is a single-cell proteomic technology using heaby metal reporter ions combined with high-dimensionsional imaging by laser ablation to mass cytometry. Using this technology, we can define cell types, spatial orientation, and functional activity of immune cell populations within colon cancers from young and late-onset patients. 

# Analysis (IN PROGRESS)
IMC images were acquired using the Hyperion system from Standard BioTools (MCD images). MCD images were visually assessed using MCD Viewer software (Standard BioTools).
Image pre-processing and segmentation were facilitated using the Steinbock toolkit (DeepCell) as previously described here: https://bodenmillergroup.github.io/steinbock/latest/cli/preprocessing/

Segmentation quality was assessed using Napari in a python script.
Downstream image and single-cell analysis was conducted in R using BioConductor packages imcRtools and Cytomapper.





