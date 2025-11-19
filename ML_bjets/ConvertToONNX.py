import os
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

import argparse
import tensorflow as tf
import tf2onnx

# tf.debugging.set_log_device_placement(True)
tf.config.set_visible_devices([], 'GPU')

def convert_h5_to_onnx(input_model_path, output_onnx_path, opset=13):
    # Load the TensorFlow model
    model = tf.keras.models.load_model(input_model_path)
    

    # Build input signature from model inputs
    input_signature = []
    for inp in model.inputs:
        input_signature.append(tf.TensorSpec(shape=inp.shape, dtype=tf.float32, name=inp.name.split(':')[0]))

    # Print input and output info
    print("Inputs:")
    for i, inp in enumerate(model.inputs):
        print(f" Input {i}: name={inp.name}, shape={inp.shape}")

    print("\nOutputs:")
    for i, out in enumerate(model.outputs):
        print(f" Output {i}: name={out.name}, shape={out.shape}")

    # Convert to ONNX
    onnx_model, _ = tf2onnx.convert.from_keras(model, input_signature=input_signature, opset=opset)
    
    # Save to output file
    with open(output_onnx_path, "wb") as f:
        f.write(onnx_model.SerializeToString())
    print(f"Model successfully converted to ONNX and saved at: {output_onnx_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a Keras .h5 model to ONNX format.")
    parser.add_argument("input_model", help="Path to the input .h5 Keras model")
    parser.add_argument("output_onnx", help="Path to save the output ONNX model")
    parser.add_argument("--opset", type=int, default=13, help="ONNX opset version (default: 13)")
    
    args = parser.parse_args()
    convert_h5_to_onnx(args.input_model, args.output_onnx, args.opset)
