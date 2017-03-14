/**@brief      Index tools for UVLM description.
 * @author     Rob Simpson
 * @contact    r.simpson11@imperial.ac.uk
 * @version    0.0
 * @date       13/05/2015
 * @pre        None
 * @warning    None
 */

#include <stdexcept>

unsigned int q_k(const unsigned int k_,
				   const unsigned int N,
				   const unsigned int no) {
	/**@brief Get lattice index of corner point based on panel index.
	 * @param k_ Panel index, starts from 0.
	 * @param N Spanwise panels.
	 * @param no Corner number within panel.
	 * @notes The equations here are from Simpson (2015) (thesis) where the
	 * indexing starts from 1, however in the implementation we need indices
	 * starting from 0, hence k = k_ + 1.
	 */

	// init k, q_k
	int k = k_ + 1;
	int q_k = 0;

	// get lattice index.
	switch(no){
	case 1 :
		q_k = (N+1)*((k-1)/N) + (k-1)%N + 1;
		break;
	case 2 :
		q_k = (N+1)*((k-1)/N) + (k-1)%N + 2;
		break;
	case 3 :
		q_k = (N+1)*((k-1)/N + 1) + (k-1)%N + 2;
		break;
	case 4 :
		q_k = (N+1)*((k-1)/N + 1) + (k-1)%N + 1;
		break;
	default :
		throw std::invalid_argument("Corner number invalid.");
	}

	return q_k-1;
}
