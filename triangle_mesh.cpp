#include <tess/triangle_mesh.h>

namespace tess
{
	// ---------------------------------------------------------------------------------------------------------------------------------------------------------
	// vertex
	// ---------------------------------------------------------------------------------------------------------------------------------------------------------

	vertex::vertex()
		: vertex(vec3::ZERO, vec3::ZERO)
	{
	}

	vertex::vertex(const vec3& _position, const vec3& _normal)
		: position(_position), normal(_normal)
	{
	}

	bool vertex::operator==(const vertex& other) const
	{
		return position == other.position && normal == other.normal;
	}

	// ---------------------------------------------------------------------------------------------------------------------------------------------------------
	// triangle_mesh
	// ---------------------------------------------------------------------------------------------------------------------------------------------------------

	bool triangle_mesh::is_valid() const
	{
		// Check if there are at least 3 vertices and 3 elements to define a triangle
		if(vertices.size() < 3 || elements.size() < 3)
		{
			return false;
		}

		// Check if mesh has at least one valid triangle
		for(unsigned int i = 0; i < elements.size() - 3; i += 3)
		{
			if(elements[i] != elements[i+1] && elements[i+1] != elements[i+2])
			{
				if(elements[i] >= vertices.size() || elements[i+1] >= vertices.size() || elements[i+2] >= vertices.size())
				{
					continue;
				}

				const auto& v1 = vertices[elements[i]].position;
				const auto& v2 = vertices[elements[i+1]].position;
				const auto& v3 = vertices[elements[i+2]].position;

				if(!math::is_valid(v1.x) || !math::is_valid(v1.y) || !math::is_valid(v1.z) ||
				   !math::is_valid(v2.x) || !math::is_valid(v2.y) || !math::is_valid(v2.z) ||
				   !math::is_valid(v3.x) || !math::is_valid(v3.y) || !math::is_valid(v3.z))
				{
					continue;
				}

				if(v1 == v2 || v2 == v3 || v1 == v3)
				{
					continue;
				}

				return true;
			}
		}

		return false;
	}
} // namespace tess
