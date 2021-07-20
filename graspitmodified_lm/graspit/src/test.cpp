#include <iostream>
#include <Python.h>

int main(int argc, char **argv)
{
    Py_Initialize();
    for (int i = 0; i < 10000; i++)
    {
        PyObject *pModule = NULL;
        PyObject *pFunc = NULL;
        PyObject *pReturn = NULL;
        PyObject *pArgs = NULL;

        std::cout << "ok1" << std::endl;
        PyRun_SimpleString("import sys");
        std::cout << "ok2" << std::endl;
        PyRun_SimpleString("sys.path.append('/root/WorkSpace/EBM_Hand')");
        std::cout << "ok3" << std::endl;
        pModule = PyImport_ImportModule("test_add");
        if (!pModule)
        {
            std::cout << "pModle is null" << std::endl;
        }

        if (PyImport_ImportModule("ebmPythonInterface") == NULL || PyErr_Occurred())
        {
            PyErr_Print();
        }
        std::cout << "ok4" << std::endl;

        if (!pModule)
        {
            std::cout << "Error: python module is null!" << std::endl;
            return 0;
        }
        else
        {
            std::cout << "python module is successful" << std::endl;
        }

        pFunc = PyObject_GetAttrString(pModule, "add");

        if (!pFunc)
        {
            std::cout << "Error: python pFunc_dof_ebm is null!" << std::endl;
            return 0;
        }
        else
        {
            std::cout << "python pFunc_dof_ebm is successful" << std::endl;
        }

        //创建参数:
        pArgs = PyTuple_New(2); //函数调用的参数传递均是以元组的形式打包的,2表示参数个数
        float a = 1.5;
        float b = 2.1;
        PyTuple_SetItem(pArgs, 0, Py_BuildValue("d", a));
        PyTuple_SetItem(pArgs, 1, Py_BuildValue("d", b));

        //返回值
        pReturn = PyEval_CallObject(pFunc, pArgs); //调用函数
        //将返回值转换为double类型
        double result = 0.0;
        PyArg_Parse(pReturn, "d", &result); //d表示转换成double型变量
        std::cout << "a + b = : " << result << std::endl;

        Py_DECREF(pArgs);
        Py_DECREF(pModule);
        Py_DECREF(pFunc);
        Py_DECREF(pReturn);
        std::cout << "i = " <<i<<std::endl;
    }
    Py_Finalize();
    return 0;
}