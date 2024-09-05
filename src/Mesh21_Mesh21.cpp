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
     */
    Mesh21::Mesh21(const std::vector<Polygon21> &cells, const std::vector<Real> &intervals): _space{cells.size()}, _time{intervals.size() - 1} {

        // Elements.
        for(Natural j = 0; j < cells.size(); ++j)
            for(Natural k = 0; k < intervals.size() - 1; ++k) {
                std::vector<Point21> points = cells[j].points();

                for(auto &point: points)
                    point += intervals[k] * 1.0_t;

                this->_elements.emplace_back(Element21{Polygon21{points}, intervals[k + 1] - intervals[k]});
            }

        // Neighbours.
        for(Natural j = 0; j < this->_space; ++j) {
            for(Natural k = 0; k < this->_time; ++k) {

                // Current element.
                Element21 current = this->_elements[j + k * this->_space];
                std::vector<Edge21> current_edges = current.b_edges();

                // Neighbours parameters.
                Integer top = (k == this->_time - 1) ? -1 :  j + (k + 1) * this->_space;
                Integer bottom = (k == 0) ? -1 : j + (k - 1) * this->_space;
                std::vector<std::array<Integer, 2>> facing(current_edges.size(), std::array<Integer, 2>{-1, -1});

                for(Natural e = 0; e < current_edges.size(); ++e) {

                    // Edge flag.
                    bool found_edge = false;

                    for(Natural i = 0; i < this->_space; ++i) {
                        for(Natural h = 0; h < this->_time; ++h) {
                            if((i == j) && (h == k))
                                continue;

                            // Candidate.
                            Element21 candidate = this->_elements[i + h * this->_space];
                            std::vector<Edge21> candidate_edges = candidate.b_edges();

                            for(Natural ce = 0; ce < candidate_edges.size(); ++ce) {
                                if(current_edges[e] == candidate_edges[ce]) {
                                    facing[e] = std::array<Integer, 2>{static_cast<Integer>(ce), static_cast<Integer>(ce)};
                                    found_edge = true;
                                    break;
                                }
                            }

                            if(found_edge)
                                break;
                        }

                        if(found_edge)
                            break;
                    }
                }
                
                this->_neighbours.emplace_back(top, bottom, facing);
            }
        }
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
     * @param j Element's space index.
     * @param k Element's time index.
     * @return std::vector<Natural> 
     */
    std::vector<Natural> Mesh21::dofs(const Natural &j, const Natural &k) const {
        #ifndef NDEBUG
        assert(j < this->_space);
        assert(k < this->_time);
        #endif

        Natural start = 0;
        std::vector<Natural> dofs;

        for(Natural h = 0; h < j * this->_time + k; ++h)
            start += this->_elements[h].dofs();

        for(Natural h = start; h < start + this->element(j, k).dofs(); ++h)
            dofs.emplace_back(h);

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