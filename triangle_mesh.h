#pragma once
#include <bl/bl.h>

namespace tess
{
	struct vertex
	{
		vertex();
		vertex(const vec3& _position, const vec3& _normal);
		bool operator==(const vertex& other) const;

		vec3 position;
		vec3 normal;
	};

	typedef unsigned int element;

	struct triangle_mesh
	{
		bool is_valid() const;
		void append(const triangle_mesh& other)
		{
			auto vertexStart = vertices.size();
			auto elementStart = elements.size();
			vertices.insert(vertices.end(), other.vertices.begin(), other.vertices.end());
			elements.insert(elements.end(), other.elements.begin(), other.elements.end());
			for(auto i = elementStart; i < elements.size(); ++i)
			{
				elements.at(i) += vertexStart;
			}
		}

		vector<vertex> vertices;
		vector<element> elements;
	};
} // namespace tess
