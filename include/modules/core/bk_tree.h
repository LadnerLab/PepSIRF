/*
 * BK-tree implementation in C++
 * Copyright (C) 2012 Eiichi Sato
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _BK_TREE_HPP_
#define _BK_TREE_HPP_

#include <map>
#include <cmath>
#include <vector>
#include <iostream>


namespace detail {

template <typename KeyType, typename MetricType, typename Distance>
class tree_node
{
private:
	typedef tree_node<KeyType, MetricType, Distance> NodeType;
	
private:
	KeyType value;
	std::map<MetricType, NodeType *> *children;

public:
	tree_node(const KeyType &key)
		: value(key), children(NULL) { }

	~tree_node() {
		if (children) {
			for (auto iter = children->begin(); iter != children->end(); ++iter)
				delete iter->second;
			delete children;
		}
	}

public:
	bool insert(NodeType *node) {
		if (!node)
			return false;

		Distance d;
		MetricType distance = d(node->value, this->value);
		if (distance == 0)
			return false; /* value already exists */

		if (!children)
			children = new std::map<MetricType, NodeType *>();

		auto iterator = children->find(distance);
		if (iterator == children->end()) {
			children->insert(std::make_pair(distance, node));
			return true;
		}

		return iterator->second->insert(node);
	}

protected:
	bool has_children() const {
		return this->children && this->children->size();
	}

protected:
	void _find_within(std::vector<std::pair<KeyType, MetricType>> &result, const KeyType &key, MetricType d) const {
		Distance f;
		MetricType n = f(key, this->value);
		if (n <= d)
			result.push_back(std::make_pair(this->value, n));

		if (!this->has_children())
			return;

		for (auto iter = children->begin(); iter != children->end(); ++iter) {
			MetricType distance = iter->first;
			if (n - d <= distance && distance <= n + d)
				iter->second->_find_within(result, key, d);
		}
	}

public:
	std::vector<std::pair<KeyType, MetricType>> find_within(const KeyType &key, MetricType d) const {
		std::vector<std::pair<KeyType, MetricType>> result;
		_find_within(result, key, d);
		return result;
	}

public:
	void dump_tree(int depth = 0) {
		for (int i = 0; i < depth; ++i)
			std::cout << "    ";
		std::cout << this->value << std::endl;
		if (this->has_children())
			for (auto iter = children->begin(); iter != children->end(); ++iter)
				iter->second->dump_tree(depth + 1);
	}
};

template <
	typename KeyType,
	typename MetricType
>
struct default_distance
{
	MetricType operator()(const KeyType &ki, const KeyType &kj) {
		return sqrt((ki - kj) * (ki - kj));
	}
};

} /* namespace detail */

template <
	typename KeyType,
	typename MetricType = double,
	typename Distance = detail::default_distance<KeyType, MetricType>
>
class bktree
{
private:
	typedef detail::tree_node<KeyType, MetricType, Distance> NodeType;

private:
	NodeType *m_top;
	size_t m_n_nodes;

public:
	bktree() : m_top(NULL), m_n_nodes(0) { }

public:
	void insert(const KeyType &key) {
		NodeType *node = new NodeType(key);
		if (!m_top) {
			m_top = node;
			m_n_nodes = 1;
			return;
		}
		if (m_top->insert(node))
			++m_n_nodes;
	};

public:
	std::vector<std::pair<KeyType, MetricType>> find_within(KeyType key, MetricType d) const {
		return m_top->find_within(key, d);
	}

	void dump_tree() {
		m_top->dump_tree();
	}

public:
	size_t size() const {
		return m_n_nodes;
	}
};

namespace bk
{
    struct hamming_distance
    {
        int operator()( const std::string& s1, const std::string& s2 )
        {
            int return_val = 0;
            std::size_t index = 0;
            for( index = 0; index < s1.length(); ++index )
                {
                    return_val += s1[ index ] != s2[ index ];
                }
            return return_val;
        }
    };
} // namespace bk

#endif /* _BK_TREE_HPP_ */
