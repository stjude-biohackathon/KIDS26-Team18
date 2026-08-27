# Biohackathon Project Template

This repository is a starting point for a three-day team project. This repository is populated with a starting template for team organization and planning. Use it to plan, build, and document work. Please adjust this repository to suit the needs of your team.

> **Team leads:** Start with the [team lead checklist](project-management/CHECKLIST.md) before the event or during your first team meeting.

## Project Profile

- **Project name:** Cell Identity in Spatial OMICs Imaging Platforms: Cluster-Based vs Cell-Based Annotation
- **Question, problem, or opportunity:** Exploring how cluster-based and cell-based annotation approaches relate to and complement each other in imaging-based spatial transcriptomics, and when each provides reliable biological insight
- **Data, inputs, or evidence:** We will be using public datasets either from GEO DataSets or from the platform companies (e.g., 10x Genomics, Bruker). Each dataset will represent either a platform or an organ type. The primary input will be an AnnData object (.h5ad file). We can use SpatialData objects (.zarr) for data visualization depending upon availability.
- **Expected output:** The expected output will be a Jupyter notebook with a straightforward analysis demonstrating the relevance of cell-based annotation. All utility functions should be compiled at a higher level so other members of the team can use the core functionalities developed during this biohackathon.
- **Tools and stack:** Python: Scanpy, AnnData, SpatialData, and scverse tools in general. Previous experience with single-cell and/or spatial transcriptomics analysis is recommended. Experience working with public datasets and repositories, such as GEO DataSets, would be helpful. Familiarity with specific organ types or biological domains (e.g., kidney diseases, cancer biology) would be a plus. Basic GitHub knowledge is also expected. None of these are mandatory.
- **Team lead:** Maycon Marção [@Mmaycon]
- **Team members and roles:** [Link to `project-management/team.md`]
- **Communication:** Slack (StJude workspace)


## Vision and Mission

- **Vision:** Explore whether cell types can be reliably identified directly from spatial transcript detections and marker genes, without relying on conventional analysis pipelines.
- **Mission:** Bring people together to test this idea across public datasets, share analytical perspectives, and introduce newcomers to spatial transcriptomics.

## About

Cell type annotation in imaging-based transcriptomics sits at the core of how we interpret biology, yet we still rely on two fundamentally different strategies without fully understanding how they relate to each other. Cluster-based annotation assigns identity at the population level, while cell-based annotation applies explicit gene rules to classify individual cells. Imaging-based technologies make this especially interesting because transcripts are directly detected as spatial spots, allowing gene expression to behave more like a robust binary signal: a transcript is either observed or not. This makes rule-based annotation intuitive and interpretable, but it also raises questions about thresholds, marker specificity, and hidden assumptions. Both approaches are widely used and often treated as interchangeable, despite little evidence that they capture the same biological reality. The goal is not to declare one method better than the other, but to understand how complementary they can be and when each provides the most reliable biological insight.

## Roadmap and Milestones

| When | Focus | Expected outcome |
| --- | --- | --- |
| Day 1 | Agree on the question, inputs, stack, roles, and first tasks | A shared plan and a first small change in the repository |
| Day 2 | Build, test, and compare approaches | A working result or clear evidence about what does not work |
| Day 3 | Stabilize, document, and present | A demo or handoff with methods, limitations, and next steps |

The goal is not a perfect production system. The goal is a clear, honest, useful result that the team can explain and others can build on.


