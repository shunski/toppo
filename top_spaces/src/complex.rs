use util::iter::Counter;
use crate::Space;
use std::collections::HashMap;
use std::collections::HashSet;
use std::cmp;

#[allow(unused)]

#[derive(Hash, Clone)] 
pub struct Simplex <LabelT> 
    where LabelT: std::clone::Clone + std::cmp::Eq + std::hash::Hash
{
    dim: usize,
    vertices_labels: Vec<LabelT>,
}

#[allow(unused)]
impl<LabelT> Simplex<LabelT> 
    where LabelT: std::clone::Clone + std::cmp::Eq + std::hash::Hash
{
    pub fn new_from(labels: Vec<LabelT>) -> Simplex<LabelT> {
        Simplex{
            dim: labels.len() - 1,
            vertices_labels: Vec::from_iter(labels.into_iter()),
        }
    }
}

pub struct SimplicialCplx <LabelT> 
    where LabelT: std::clone::Clone + std::cmp::Eq + std::hash::Hash
{
    num_vertices: usize,
    maximal_cells: Vec<Vec<usize>>,
    vertices_labeling: HashMap<LabelT, usize>,
    connected_components: Vec<Vec<usize>>
}

#[allow(unused)]
impl<LabelT> SimplicialCplx<LabelT> 
    where LabelT: std::clone::Clone + std::cmp::Eq + std::hash::Hash
{
    pub fn new() -> SimplicialCplx<LabelT> {
        SimplicialCplx{
            num_vertices: 0,
            maximal_cells: Vec::new(),
            vertices_labeling: HashMap::new(),
            connected_components: Vec::new(),
        }
    }

    pub fn new_from_zero_skeleton(vertices_labels: Vec<LabelT>) -> SimplicialCplx<LabelT> {
        SimplicialCplx {
            num_vertices: vertices_labels.len(),
            maximal_cells: 
                Counter::new_with_max(vertices_labels.len())
                .map(|x| vec![x])
                .collect(),
            connected_components: 
                Counter::new_with_max(vertices_labels.len())
                .map(|x| vec![x])
                .collect(),
            vertices_labeling: 
                vertices_labels
                .into_iter()
                .zip(Counter::new())
                .collect(),
        }
    }

    pub fn new_from(simplex: Simplex<LabelT>) -> SimplicialCplx<LabelT> 
    {
        let mut cplx = SimplicialCplx::new();
        cplx.attach(simplex);
        cplx
    }

    fn get_connected_component_id(&self, vertex: usize) ->  usize {
        self.connected_components
            .iter().enumerate()
            .find(|(_, component)| component.contains(&vertex))
            .unwrap()
            .0
    }

    pub fn attach(&mut self, simplex: Simplex<LabelT>){
        // add vertices to the complex if necessary
        for v_label in simplex.vertices_labels.iter() {
            if !self.vertices_labeling.contains_key(&v_label) {
                self.vertices_labeling.insert(v_label.clone(), self.num_vertices);
                self.maximal_cells.push(vec![self.num_vertices]);
                self.connected_components.push(vec![self.num_vertices]);
                self.num_vertices += 1;
            }
        }

        // build the ordered simplex using Simplex data
        let mut simplex:Vec<usize> = simplex
            .vertices_labels
            .into_iter()
            .map(|x| match self.vertices_labeling.get(&x) {
                Some(&vertex) => vertex, 
                None => panic!(), 
            })
            .collect();
        simplex.sort();

        // if the new simplex is contained in the complex then do nothing
        if self.maximal_cells.iter().any(|max_cell| simplex.iter().all(|&vertex| max_cell.contains(&vertex))) {
            return;
        }

        // Otherwise we will add this simplex to the complex
        // First remove all the simplices that are contained in the new simplex, 
        // then modify the connected components (some components may be glued by adding),
        // and finally add the new simplex,.
        self.maximal_cells.retain(|cell| cell.iter().all(|vertex| simplex.contains(&vertex)));

        let mut components_to_be_removed = 
            simplex.iter()
            .skip(1)
            .map(|&vertex| self.get_connected_component_id(vertex))
            .collect::<Vec<_>>();
        components_to_be_removed.sort_by(|a, b| b.cmp(a));
        components_to_be_removed.dedup();
        let component_to_be_merged = self.get_connected_component_id(simplex[0]);
        Counter::new_with_max(self.num_vertices).for_each(
            |v_idx| {
                if components_to_be_removed.iter().any(|&component_idx| self.connected_components[component_idx].contains(&v_idx)) {
                    self.connected_components[component_to_be_merged].push(v_idx);
                }
            }
        );
        self.connected_components[component_to_be_merged].sort();
        self.connected_components[component_to_be_merged].dedup();
        components_to_be_removed.iter().for_each(|&component| {self.connected_components.remove(component); });

        self.maximal_cells.push(simplex);
    }

    pub fn get_labeled_zero_skeleton(&self) -> HashSet<LabelT> {
        return self.vertices_labeling.keys().cloned().collect();
    }

    pub fn get_dim(&self) -> usize {
        self.maximal_cells
            .iter()
            .fold(0, |current_max, cell| cmp::max(current_max, cell.len()-1))
    }

    pub fn is_connected(&self) -> bool {
        return self.connected_components.len() == 1
    }
}




impl<LabelT> Space for SimplicialCplx<LabelT> 
    where LabelT: std::clone::Clone + std::cmp::Eq + std::hash::Hash
{
    fn compute_homology(&self) -> Vec<usize> {
        
    }
}


#[cfg(test)]
mod simplicial_cplx_tests {
    #[test]
    fn object_instantiation() {
        use crate::copmlex::*;

        let mut cplx = SimplicialCplx::new_from_zero_skeleton(vec!["v0".to_string(), "v1".to_string()]);
        assert_eq!(cplx.is_connected(), false);

        let mut labels = vec!["v1".to_string(), "v2".to_string(), "v3".to_string()];
        let simplex = Simplex::new_from(labels.clone());
        assert_eq!(simplex.dim, 2);

        cplx.attach(simplex.clone());
        labels.push("v0".to_string());
        assert_eq!(
            cplx.get_labeled_zero_skeleton(), 
            HashSet::from_iter(labels.clone().into_iter())
        );

        cplx.attach(simplex.clone());
        assert_eq!(
            cplx.get_labeled_zero_skeleton(), 
            HashSet::from_iter(labels.clone().into_iter())
        );

        let mut labels = vec!["v0".to_string(), "v1".to_string()];
        let simplex = Simplex::new_from(labels.clone());
        let cplx = SimplicialCplx::new_from(simplex);
        assert_eq!(cplx.get_dim(), 1);

        labels.pop();
        labels.push("v2".to_string());
        labels.push("v3".to_string());
        labels.push("v4".to_string());
        let simplex = Simplex::new_from(labels.clone());
        let cplx = SimplicialCplx::new_from(simplex);
        assert_eq!(cplx.get_dim(), 3);
        assert_eq!(cplx.is_connected(), true);
    }
}