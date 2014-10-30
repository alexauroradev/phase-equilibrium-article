#ifndef __UTIL_H__
#define __UTIL_H__

#include <string>
#include <cstdlib>
#include <cstdio>
#include <malloc.h>
#include <cstdarg>

#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>
#include <map>

#include <vector>
#include <cstddef>

#include <stdexcept>

namespace util {

static inline std::string format(const char *fmt, ...)
#if defined(__GNUC__) || defined(__GNUG__) || defined (__clang__)
	__attribute__((format (printf, 1, 2)));
#else
	;
#endif

static inline std::string format(const char *fmt, ...)
{
	size_t size = 256;
	char *buf = (char *)malloc(size);
	va_list args;

	va_start(args, fmt);
	int ret = vsnprintf(buf, size, fmt, args);
	va_end(args);

	std::string tmp;

	if (ret < 0)
		tmp = "";
	else if ((size_t)ret >= size) {
		size = ret + 1;
		buf = (char *)realloc(buf, size);
		if (buf) {
			va_start(args, fmt);
			ret = vsnprintf(buf, size, fmt, args);
			va_end(args);
			if (ret < 0)
				tmp = "";
			else
				tmp = buf;
		} else
			tmp = "";
	} else
		tmp = buf;

	free(buf);
	return std::string(tmp);
}

// trim from start
static inline std::string &ltrim(std::string &s) {
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(isspace))));
	return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
	s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(isspace))).base(), s.end());
	return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
	return ltrim(rtrim(s));
}

template<typename Key, typename Value>
const Value &get(const std::map<Key, Value> &map, const Key &key) {
	typename std::map<Key, Value>::const_iterator pos = map.find(key);
	if (pos == map.end())
		throw std::out_of_range("Invalid key value for const map access");
	return pos->second;
}

template<class T>
class ptr_vector {
	std::vector<T *> m;
	public:
		explicit ptr_vector(const size_t n) : m(n, 0) { }

	size_t size() const {
		return m.size();
	}

	size_t capacity() const {
		return m.capacity();
	}

	T *put(const ptrdiff_t i, T *v) {
		T *old = m[i];
		m[i] = v;
		return old;
	}

	bool is_null(const ptrdiff_t i) {
		return m[i] == 0;
	}

	T &operator[](const ptrdiff_t i) {
		return *m[i];
	}

	const T &operator[](const ptrdiff_t i) const {
		return *m[i];
	}

	~ptr_vector() {
		for (typename std::vector<T *>::const_iterator i = m.begin(); i != m.end(); ++i)
			delete *i;
	}
};

}

#endif
