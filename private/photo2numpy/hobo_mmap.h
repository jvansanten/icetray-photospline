/*
 * While Python has an mmap interface, it's difficult work with from C.
 * POSIX is all the portability we need.
 */

#include <Python.h>
#include <sys/mman.h>
#include <errno.h>

typedef struct {
	PyObject_HEAD
	void *pages;
	size_t length;
} mmap_wrapper;

static void
mmap_wrapper_dealloc(mmap_wrapper *self)
{
	munmap(self->pages, self->length);	
}

static PyTypeObject mmap_wrapper_type = {
	PyObject_HEAD_INIT(NULL)
	0,                         /*ob_size*/
	"mmap_wrapper",             /*tp_name*/
	sizeof(mmap_wrapper), /*tp_basicsize*/
	0,                         /*tp_itemsize*/
	(destructor)mmap_wrapper_dealloc, /*tp_dealloc*/
	0,                         /*tp_print*/
	0,                         /*tp_getattr*/
	0,                         /*tp_setattr*/
	0,                         /*tp_compare*/
	0,                         /*tp_repr*/
	0,                         /*tp_as_number*/
	0,                         /*tp_as_sequence*/
	0,                         /*tp_as_mapping*/
	0,                         /*tp_hash */
	0,                         /*tp_call*/
	0,                         /*tp_str*/
	0,                         /*tp_getattro*/
	0,                         /*tp_setattro*/
	0,                         /*tp_as_buffer*/
	Py_TPFLAGS_DEFAULT,        /*tp_flags*/
	"Thin wrapper around mmap'd pages",           /* tp_doc */
};

static mmap_wrapper*
mmap_wrapper_new(int fd, off_t offset, size_t length, void **start)
{
	void *pages;
	off_t front_padding, map_start;
	size_t pagesize, map_length;
	mmap_wrapper *self;

	/* Align the requested mapping to page boundaries. */
	pagesize = getpagesize();
	map_start = pagesize * (offset / pagesize);
	front_padding = offset - map_start;
	map_length = front_padding + length;
	if (map_length % pagesize != 0) {
		size_t min_length = map_length;
		map_length += pagesize - (min_length % pagesize);
	}

	/* Map it! */
	pages = mmap(NULL, map_length, PROT_READ | PROT_WRITE,
	    MAP_FILE | MAP_PRIVATE, fd, map_start);

	if (pages == MAP_FAILED) {
		PyErr_SetFromErrno(PyExc_RuntimeError);
		return (NULL);
	}

	/* Where does the requested data actually start? */
	*start = (char*)pages + front_padding;

	self = (mmap_wrapper*)mmap_wrapper_type.tp_alloc(&mmap_wrapper_type, 0);
	self->pages = pages;
	self->length = map_length;

	return (self);
}

