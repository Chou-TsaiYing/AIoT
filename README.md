# Federated learning with intelligent edge server in  massive IoT scenario

## Introduction

We consider three diffenet resource allocation methods in No shared (Small-private) service scenarios. Three different resource allocation methods are **Average allocation** , **Cost bisection**, and **Optimize**. Under average training delay constraint and average inferencing delay constraint, we aim to minimize the total weighted computation capacity. Finally, we will compare the three resource allocate methods.

## Build With
* C/C++

## Program Architecture and Description
### Architecure


| Program File         | Function                                                                                                                                                                                                                     |
| -------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
|       twister.c               |    generates a random number on (0,1)-real-interval                                                                                                                                                                                                                          |
| ML_ResourceAlloc.cpp | Three diffenet resource allocation methods in No shared (Small-private) service scenarios.Under average training delay constraint and average inferencing delay, we aim to minimize the total weighted computation capacity. |

### Description

#### **How to run**
```sh
    g++ -o .\ML_ResourceAlloc .\ML_ResourceAlloc.cpp
```
```sh
    .\ML_ResourceAlloc
```

Enter the number to choose one of algorithm

![](https://i.imgur.com/aXopP3S.png)

#### **Output File**

Method: **Average allocation**

| File                 | Content                                                                                                                                       |
| -------------------- | --------------------------------------------------------------------------------------------------------------------------------------------- |
| base_capacity.csv    | Computation node capacity, total capacity, total weighted capacity during training and prediction in 18 cases of No Shared service scenarios with Average allocation Method. |
| average_C_node_P.csv | Allocated capacity for each node with Average allocation Method during prediction.                                                                     |
| average_C_node_T.csv | Allocated capacity for each node with Average allocation Method during training.                                                                     |

Method: **Cost bisection**

| File                 | Content                                                                                                                                       |
| -------------------- | --------------------------------------------------------------------------------------------------------------------------------------------- |
| bisection_capacity.csv    | Computation node capacity, total capacity, total weighted capacity during training and prediction in 18 cases of No Shared service scenarios with Cost bisection method. |
| bisection_C_node_P.csv | Allocated capacity for each node  with Cost bisection method during prediction.                                                                     |
| bisection_C_node_T.csv | Allocated capacity for each node with Cost bisection method during training.                                                                     |

Method: **Optimize**

| File                 | Content                                                                                                                                       |
| -------------------- | --------------------------------------------------------------------------------------------------------------------------------------------- |
| optimizeAlloc_capacity.csv   | Computation node capacity, total capacity, total weighted capacity during training and prediction in 18 cases of No Shared service scenarios with Optimize method. |
| bisection_C_node_P.csv | Allocated capacity for each node with Optimize method during prediction                                                                     |
| bisection_C_node_T.csv | Allocated capacity for each node with Optimize method during training 

#### **Result**
Total weighted capacity compare

Optimize method is better than Average allocation method and Cost bisection method.
![](https://i.imgur.com/yL7fDZd.png)
