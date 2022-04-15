#Self-supervised learning for analysis of temporal and morphological drug effects in cancer cell imaging data

## Abstract

In this work, we propose two novel methodologies to study temporal and morphological phenotypic effects caused by different experimental conditions using imaging data. As a proof of concept, we apply them to analyze drug effects in 2D cancer cell cultures. We train a convolutional autoencoder on 1M images dataset with random augmentations and multi-crops to use as feature extractor. We systematically compare it to the pretrained state-of-the-art models. We further use the feature extractor in two ways. First, we apply distance-based analysis and dynamic time warping to cluster temporal patterns of 31 drugs. We identify clusters allowing annotation of drugs as having cytotoxic, cytostatic, mixed or no effect. Second, we implement an adversarial/regularized learning setup to improve classification of 31 drugs and visualize image regions that contribute to the improvement. We increase top-3 classification accuracy by 8% on average and mine examples of morphological feature importance maps. We provide the feature extractor and the weights to foster transfer learning applications in biology. We also discuss utility of other pretrained models and applicability of our methods to other types of biomedical data.

<img src="https://github.com/mauromiguelm/toxicity-classifier/blob/master/img/distance_methods.png" alt="Methods developed in the paper." title="Methods developed in the paper."
width="600"/>

<img src="https://github.com/mauromiguelm/toxicity-classifier/blob/master/img/MIDL_lastfigure_halfcropped.png" alt="Outcomes of time-series analysis of encodings." title="Outcomes of time-series analysis of encodings."
width="750"/>

##About
For more details, please check the manuscript:
Dmitrenko, A., Masiero, M. M., & Zamboni, N. (2022). Self-supervised learning for analysis of temporal and morphological drug effects in cancer cell imaging data. Proceedings of Machine Learning Research-Under Review, 1–17. https://doi.org/10.48550/arxiv.2203.04289
