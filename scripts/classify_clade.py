import re
import csv
from collections import OrderedDict
import argparse

classification = {}   ##note the result of classification

class TreeNode:
    def __init__(self, name=None):
        self.name = name
        self.children = []  # 用于存储子节点
        self.value = None

    def add_child(self, child_node):
        self.children.append(child_node)

def read_fasta(filename):
    with open(filename, 'r') as file:
        content = file.read()
        return content

def is_english_name(name):
    return bool(re.match(r'^[A-Za-z]+[0-9]*$', name))

def parse_newick(data):
    stack = []
    current_node = TreeNode()
    root = current_node
    i = 0

    # 去除最外层括号
    if data.startswith('(') and data.endswith(')'):
        data = data[1:-1]

    node_pattern = re.compile(r"('(\w+)':([\d.e-]+))")
    edge_pattern = re.compile(r"(\d+\.\d+):(\d+\.\d+)")

    while i < len(data):
        if data[i] == '('  or data[i:i+2] == ',(':
            if data[i:i+2] == ',(':
                i+=1
            # 创建新节点并进入下一层
            new_node = TreeNode()
            
            stack.append(current_node)
            
            current_node.add_child(new_node)
            current_node = new_node
            i += 1
        if data[i:i+2] == ',\'':
            i+=1
        
            new_node = TreeNode()
            current_node.add_child(new_node)
            match = node_pattern.match(data, i)
            if match:
                node_name = match.group(2)
                node_value = float(match.group(3))
                new_node.value = node_value
                new_node.name = node_name
                i = match.end()
            
        elif data[i] == ')':
            # 返回上一层
            current_node = stack.pop()
            match_tmp = re.match(r'\)\d*\.?\d*:\d*\.?\d*', data[i:])
            if match_tmp:
                i += len(match_tmp.group(0))
                if i+1<len(data) and data[i+1] == ';':
                    i+=1
        else:
            # read the value of node
            match = node_pattern.match(data, i)
            if match:
                tmp_node = TreeNode()
                node_name = match.group(2)
                node_value = float(match.group(3))
                tmp_node.value = node_value
                tmp_node.name = node_name
                
                current_node.add_child(tmp_node)
                
                i = match.end()
            if data[i] == ';':
                break
        

    return root


def collect_label_paths(node, path=[], label_paths={}):
    if node is None:
        return label_paths

    # 更新当前路径
    current_path = path + [node]

    # 如果是标签节点，记录路径
    if node.name is not None and is_english_name(node.name):
            label_paths[node.name] = current_path

    # 继续遍历子节点
    for child in node.children:
        collect_label_paths(child, current_path, label_paths)

    return label_paths

def classify_node(node, label_paths):
    if node is not None and is_english_name(node[-1].name):
        return None  # 忽略标签节点或空节点

    # 计算与每个叶子节点标签路径的最大重叠
    max_overlap = 0
    best_label = None
    for label, path in label_paths.items():
        overlap = sum(1 for i, j in zip(node, path) if i == j)
        if overlap > max_overlap:
            max_overlap = overlap
            best_label = label

    return best_label


def DFS(root, label_paths, current_path):
    if root is None:
        return

    # 在当前路径中添加此节点
    current_path.append(root)

    # 如果是叶子节点，执行分类
    if not root.children and not is_english_name(root.name):
        classification[root.name] = classify_node(current_path, label_paths)
    else:
        # 递归调用子节点
        for child in root.children:
            DFS(child, label_paths, current_path.copy())

    # 回溯时移除当前节点
    current_path.pop()

def main(data_file, output_file):
    # 读取数据
    data = read_fasta(data_file)
    print("start to build tree")
    root = parse_newick(data)
    print("build tree successfully")
    # 收集标签路径
    label_paths = collect_label_paths(root)
    my_path = []
    DFS(root, label_paths, my_path)
    print("classify successfully")

    # 修改排序逻辑，将 None 视为 'Unsorted'
    sorted_classified_result = OrderedDict(sorted(
        classification.items(), 
        key=lambda x: x[1] if x[1] is not None else 'Unsorted'
    ))

    # print(sorted_classified_result)
    
    # save to file
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Digit Node', 'Class'])  # 写入表头
        for key, value in sorted_classified_result.items():
            writer.writerow([key, value if value is not None else 'Unsorted'])  # 将 None 替换为 'Unsorted'

    print(f"The result has been saved to {output_file} file.")
    
if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Process tree classification.')
    parser.add_argument('data_file', type=str, help='Input data file path (e.g., ./cin_copia.rt.tree)')
    parser.add_argument('output_file', type=str, help='Output CSV file path (e.g., classification_results.csv)')
    
    args = parser.parse_args()
    main(args.data_file, args.output_file)