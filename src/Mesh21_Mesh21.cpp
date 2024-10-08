/**
 * @file Mesh21_Mesh21.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Mesh21/Mesh21.hpp implementation.
 * @date 2024-07-25
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    // Constructors.

    /**
     * @brief Default constructor.
     * 
     * @param cells Space cells.
     * @param intervals Time intervals, length = N + 1: [t0, t1, t2, ..., tN].
     * @param p Space degree.
     * @param q Time degree.
     */
    Mesh21::Mesh21(const std::vector<Polygon21> &cells, const std::vector<Real> &intervals, const Natural &p, const Natural &q): _space{cells.size()}, _time{intervals.size() - 1} {
        
        #ifndef NVERBOSE
        std::cout << "[Ivo] Mesh21" << std::endl;
        std::cout << "\t[Mesh21] Evaluating elements" << std::endl;
        #endif

        // Elements.
        for(Natural j = 0; j < this->_time; ++j)
            for(Natural k = 0; k < this->_space; ++k) {
                std::vector<Point21> points = cells[k].points();

                for(auto &point: points)
                    point += intervals[j] * 1.0_t;

                this->_elements.emplace_back(Element21{Polygon21{points}, intervals[j + 1] - intervals[j], p, q});
            }

        #ifndef NVERBOSE
        std::cout << "\t[Mesh21] Evaluating neighbours" << std::endl;
        #endif

        // Neighbours.
        for(Natural j = 0; j < this->_time; ++j)
            for(Natural k = 0; k < this->_space; ++k) {

                // Current element.
                Element21 current = this->_elements[j * this->_space + k];
                std::vector<Edge21> current_edges = current.b_edges();

                // Neighbours parameters.
                Integer top = (j == this->_time - 1) ? -1 :  (j + 1) * this->_space + k;
                Integer bottom = (j == 0) ? -1 : (j - 1) * this->_space + k;
                std::vector<std::array<Integer, 2>> facing(current_edges.size(), std::array<Integer, 2>{-1, -1});

                for(Natural e = 0; e < current_edges.size(); ++e) {

                    // Edge flag.
                    bool found_edge = false;

                    for(Natural i = 0; i < this->_space; ++i) {
                        if(i == k)
                            continue;

                        // Candidate.
                        Element21 candidate = this->_elements[j * this->_space + i];
                        std::vector<Edge21> candidate_edges = candidate.b_edges();

                        for(Natural ce = 0; ce < candidate_edges.size(); ++ce) {
                            if(current_edges[e] == candidate_edges[ce]) {
                                facing[e][0] = static_cast<Integer>(j * this->_space + i);
                                facing[e][1] = static_cast<Integer>(ce);

                                found_edge = true;
                                break;
                            }
                        }

                        if(found_edge)
                            break;
                    }
                }
                
                this->_neighbours.emplace_back(top, bottom, facing);
            }

        #ifndef NVERBOSE
        std::cout << "\t[Mesh21] Exited" << std::endl;
        #endif
    }

    // Parameters.

    /**
     * @brief Mesh' dofs.
     * 
     * @return Natural 
     */
    Natural Mesh21::dofs() const {
        Natural dofs = 0;

        for(const auto &element: this->_elements)
            dofs += element.dofs();

        return dofs;
    }

    /**
     * @brief Global to local dofs.
     * 
     * @param j Element's index.
     * @return std::vector<Natural> 
     */
    std::vector<Natural> Mesh21::dofs(const Natural &j) const {
        #ifndef NDEBUG
        assert(j < this->_elements.size());
        #endif

        Natural start = 0;
        std::vector<Natural> dofs;

        for(Natural k = 0; k < j; ++k)
            start += this->_elements[k].dofs();

        for(Natural k = start; k < start + this->element(j).dofs(); ++k)
            dofs.emplace_back(k);

        return dofs;
    }

    /**
     * @brief Global to local dofs.
     * 
     * @param j Time slab's index.
     * @return std::vector<Natural> 
     */
    std::vector<Natural> Mesh21::dofs_t(const Natural &j) const {
        #ifndef NDEBUG
        assert(j < this->_time);
        #endif

        Natural counter = 0;
        std::vector<Natural> dofs;

        for(Natural k = 0; k < j * this->_space; ++k)
            counter += this->_elements[k].dofs();

        for(Natural k = 0; k < this->_space; ++k) {
            for(Natural h = 0; h < this->element(j * this->_space + k).dofs(); ++h) {
                dofs.emplace_back(counter);
                ++counter;
            }
        }

        return dofs;
    }

    // Output.

    /**
     * @brief Mesh output. CSV format.
     * 
     * @param ost 
     * @param mesh Mesh.
     * @return std::ostream& 
     */
    std::ostream &operator <<(std::ostream &ost, const Mesh21 &mesh) {
        for(auto it = mesh._elements.begin(); it < mesh._elements.end() - 1; ++it)
            ost << *it << std::endl;
        
        return ost << *--mesh._elements.end() << std::flush;
    }

}