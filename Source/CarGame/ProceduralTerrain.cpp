#include "ProceduralTerrain.h"
#include "DrawDebugHelpers.h"
#include "Misc/DateTime.h"
#include "Misc/Timespan.h"
#include "Containers/UnrealString.h"
#include "ProceduralMeshComponent.h"
#include "Async/Async.h"
#include "Engine/World.h"            
#include "PhysicsEngine/BodySetup.h" 
#include "PhysicsEngine/ConvexElem.h" 
#include "Containers/Array.h"       
#include "Math/Vector.h"           
#include "Math/Vector2D.h"           
#include "Math/Color.h"            
#include "Components/SceneComponent.h"
#include "GameFramework/Character.h"
#include "GameFramework/Pawn.h"
#include "GameFramework/PlayerController.h"
#include "GameFramework/Actor.h"
#include "GameFramework/Controller.h"
#include "GameFramework/MovementComponent.h"
#include "Components/InputComponent.h"
#include "GameFramework/CharacterMovementComponent.h"


AProceduralTerrain::AProceduralTerrain()
{
	PrimaryActorTick.bCanEverTick = true;
}

// Called when the game starts or when spawned
void AProceduralTerrain::BeginPlay()
{
	Super::BeginPlay();
	randomOffset.Value = rand() % 100000;
	randomOffset.Key = rand() % 100000;
	MeshComponent.Empty();
	TerrainComponents.Empty();
	PlayerPawn = GetWorld()->GetFirstPlayerController()->GetPawn();
	currentGridSize = Width * Scale;
	PlayerGridPos = FVector2D(0, 0);
	RealtimeMesh = GetRealtimeMeshComponent()->InitializeRealtimeMesh<URealtimeMeshSimple>();
	RealtimeMesh->SetupMaterialSlot(0, "PrimaryMaterial", TerrainMaterial);
	FRealtimeMeshCollisionConfiguration CollisionConfig;
	CollisionConfig.bUseAsyncCook = true;  // Allows async cooking for better performance
	CollisionConfig.bUseComplexAsSimpleCollision = true;
	CollisionConfig.bShouldFastCookMeshes = true;
	RealtimeMeshComponent->SetMobility(EComponentMobility::Static);
	RealtimeMeshComponent->SetGenerateOverlapEvents(false);
	RealtimeMeshComponent->SetCollisionProfileName(UCollisionProfile::BlockAll_ProfileName);
	RealtimeMeshComponent->SetCollisionEnabled(ECollisionEnabled::QueryAndPhysics);
	SetRootComponent(RealtimeMeshComponent);
	RealtimeMesh->SetCollisionConfig(CollisionConfig);
	GenerateTerrain();

}


void AProceduralTerrain::GenerateTerrain()
{
	FDateTime StartTime = FDateTime::Now();
	GeneratePath();
	SmoothPathPointsHeight(PathHeightSmooth);
	GeneratePathMesh();
	DisplayPathMesh();
	UpdateTerrain();

	FDateTime EndTime = FDateTime::Now();
	FTimespan Duration = EndTime - StartTime;
	UE_LOG(LogTemp, Display, TEXT("TERRAIN GEN COMPLETE IN %f SECONDS"), Duration.GetTotalSeconds());
}

// Called every frame
void AProceduralTerrain::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);
	if (GetWorld()->GetFirstPlayerController()->IsInputKeyDown(EKeys::M)) {
		ResetPlayerPhysics();
	}

	FVector PlayerPos = PlayerPawn->GetActorLocation();
	if (PlayerPos.X > ((PlayerGridPos.X * currentGridSize) + (currentGridSize))) {
		//right
		PlayerGridPos.X++;
		UpdateTerrain();
		UE_LOG(LogTemp, Display, TEXT("RIGHT"));
	}
	else if (PlayerPos.X < ((PlayerGridPos.X * currentGridSize))) {
		//left
		PlayerGridPos.X--;
		UpdateTerrain();
		UE_LOG(LogTemp, Display, TEXT("LEFT"));
	}
	else if (PlayerPos.Y > ((PlayerGridPos.Y * currentGridSize) + currentGridSize)) {
		//top
		PlayerGridPos.Y++;
		UpdateTerrain();
		UE_LOG(LogTemp, Display, TEXT("TOP"));
	}
	else if (PlayerPos.Y < ((PlayerGridPos.Y * currentGridSize))) {
		//bottom
		PlayerGridPos.Y--;
		UpdateTerrain();
		UE_LOG(LogTemp, Display, TEXT("BOTTOM"));
	}

}

void AProceduralTerrain::ResetPlayerPhysics()
{
	if (!PlayerPawn) return;

	// Reset rotation with Teleport flag
	PlayerPawn->SetActorRotation(FRotator::ZeroRotator, ETeleportType::TeleportPhysics);

	if (ACharacter* Character = Cast<ACharacter>(PlayerPawn))
	{
		if (UCharacterMovementComponent* MovementComp = Character->GetCharacterMovement())
		{
			MovementComp->Velocity = FVector::ZeroVector;
		}
	}

	if (UPrimitiveComponent* RootPrimitive = Cast<UPrimitiveComponent>(PlayerPawn->GetRootComponent()))
	{
		RootPrimitive->SetPhysicsLinearVelocity(FVector::ZeroVector);
		RootPrimitive->SetPhysicsAngularVelocityInDegrees(FVector::ZeroVector);
	}
}




void AProceduralTerrain::UpdateTerrain() {
	TArray<TerrainComponent*> terrainsToBeMoved;
	TArray<TerrainComponent*> ComponentsToUpdate;
	TArray<std::pair<FVector2D, int>> terrainsToMake; //desired grid pos and LOD
	float range4 = 16.0f;
	int distanceSearch = range4 + 4;
	for (int x = -distanceSearch; x <= distanceSearch; x++) {
		for (int y = -distanceSearch; y <= distanceSearch; y++) {
			FVector2D GridPosition(x, y);
			GridPosition += PlayerGridPos;
			float distance = FVector2D::Distance(GridPosition, PlayerGridPos);

			if (distance <= 2.0f)
			{
				// High LOD
				TerrainComponent* component = FindTerrainComponent(GridPosition);
				if (component != nullptr && component->GetLOD() != 0) {
					component->SetLOD(0);
					terrainsToBeMoved.Add(component);
					terrainsToMake.Add(std::make_pair(GridPosition, 0));
				}
				else if (component == nullptr) {
					terrainsToMake.Add(std::make_pair(GridPosition, 0));
				}
			}
			else if (distance <= 6.0f)
			{
				// Medium LOD
				TerrainComponent* component = FindTerrainComponent(GridPosition);
				if (component != nullptr && component->GetLOD() != 1) {
					component->SetLOD(1);
					terrainsToBeMoved.Add(component);
					terrainsToMake.Add(std::make_pair(GridPosition, 1));
				}
				else if (component == nullptr) {
					terrainsToMake.Add(std::make_pair(GridPosition, 1));
				}
			}
			else if (distance <= 10.0f)
			{
				// Low LOD
				TerrainComponent* component = FindTerrainComponent(GridPosition);
				if (component != nullptr && component->GetLOD() != 3) {
					component->SetLOD(3);
					terrainsToBeMoved.Add(component);
					terrainsToMake.Add(std::make_pair(GridPosition, 3));
				}
				else if (component == nullptr) {
					terrainsToMake.Add(std::make_pair(GridPosition, 3));
				}
			}
			else if (distance <= range4)
			{
				// Lowest LOD
				TerrainComponent* component = FindTerrainComponent(GridPosition);
				if (component != nullptr && component->GetLOD() != 4) {

					component->SetLOD(4);
					terrainsToBeMoved.Add(component);
					terrainsToMake.Add(std::make_pair(GridPosition, 4));
				}
				else if (component == nullptr) {
					terrainsToMake.Add(std::make_pair(GridPosition, 4));
				}
			}
			else {
				TerrainComponent* component = FindTerrainComponent(GridPosition);
				if (component != nullptr) {
					terrainsToBeMoved.Add(component);
				}

			}
		}
	}

	UE_LOG(LogTemp, Display, TEXT("Terrains being moved: %i"), terrainsToBeMoved.Num());
	UE_LOG(LogTemp, Display, TEXT("Terrains required: %i"), terrainsToMake.Num());

	int terrainPiecesMade = 0;
	for (std::pair < FVector2D, int> terrainToMake : terrainsToMake) {
		bool foundMatch = false;
		for (TerrainComponent* terrainToMove : terrainsToBeMoved) {
			if (terrainToMove->GetLOD() == terrainToMake.second) {
				foundMatch = true;
				terrainToMove->SetGridPosition(terrainToMake.first);
				ComponentsToUpdate.Add(terrainToMove);
				terrainsToBeMoved.Remove(terrainToMove);
				break;
			}
		}
		if (!foundMatch) {
			ComponentsToUpdate.Add(CreateTerrainComponent(terrainToMake.first, terrainToMake.second));
			terrainPiecesMade++;
		}

	}
	UE_LOG(LogTemp, Display, TEXT("New Terrains: %i"), terrainPiecesMade);

	ParallelFor((ComponentsToUpdate.Num()), [this, ComponentsToUpdate](int32 Index)
		{
			GenerateTerrainSection(ComponentsToUpdate[Index]);
		});

}


TerrainComponent* AProceduralTerrain::FindTerrainComponent(const FVector2D& GridPosition) {
	for (TerrainComponent* Component : TerrainComponents)
	{
		if (Component && Component->GetGridPosition() == GridPosition)
		{
			return Component;
		}
	}
	return nullptr;
}

TerrainComponent* AProceduralTerrain::CreateTerrainComponent(const FVector2D& GridPosition, int LOD)
{
	totalCreated++;
	FString ComponentName = FString::Printf(TEXT("MeshComponent%d"), TerrainComponents.Num());
	TerrainComponent* NewComponent = new TerrainComponent(nullptr, GridPosition, LOD, TerrainComponents.Num()); //nullptr
	NewComponent->SetIsActive(true);
	TerrainComponents.Add(NewComponent);
	return NewComponent;
}


float AProceduralTerrain::CalculateNoiseAtPoint(int32 X, int32 Y) const {
	float FinalHeight = 0.0f;
	for (const FTerrainLayer& Layer : TerrainLayers)
	{
		if (Layer.Enabled) {
			float LayerHeight = FMath::PerlinNoise2D(FVector2D(X + randomOffset.Value, Y + randomOffset.Key) * (Layer.Frequency / 4)) * Layer.Amplitude;
			FinalHeight += LayerHeight * Layer.Contribution;
		}
	}
	return FinalHeight;
}


float AProceduralTerrain::CalculateHeightOnPath(int32 X, int32 Y) const
{

	float closestDist = INFINITY;
	FVector closestPoint;
	int MaxSearchRadius = 1000;
	FIntPoint Cell(FMath::FloorToInt(X * Scale / GridSize), FMath::FloorToInt(Y * Scale / GridSize));
	FVector p = FVector(X * Scale, Y * Scale, 0);
	int32 SearchRadius = 1; // Start with a small radius
	bool PointFound = false;
	for (FVector point : PathVertices) {
		if (FVector::DistXY(p, point) < closestDist) {
			closestDist = FVector::Dist2D(p, point);
			closestPoint = point;
		}
	}
	return closestPoint.Z - (HeightAdjust * 2);
}


float AProceduralTerrain::CalculateHeight(int32 X, int32 Y) const
{
	float pointHeight = CalculateNoiseAtPoint(X, Y);
	if (IsOnPath(X, Y, true))
	{
		return CalculateHeightOnPath(X, Y);
	}
	else if (IsOnPath(X, Y, false)) {
		//return AverageHeightNeighbours(X, Y, SmoothingSize);
		float pathHeight = CalculateHeightOnPath(X, Y);
		float lerpAmount = (DistFromPath(X, Y, false) / (SmoothingThicknessOffset - ThicknessOffset));
		return FMath::Lerp(pathHeight, pointHeight, lerpAmount);
	}
	return pointHeight;
}


bool AProceduralTerrain::IsOnPath(int32 X, int32 Y, bool useOffset) const
{
	FVector2D Point2D(X * Scale, Y * Scale);
	FVector Point(Point2D.X, Point2D.Y, 0.0f);
	for (int32 i = 0; i < PathPoints.Num(); ++i)
	{
		// Get current and next point, converting to FVector
		FVector2D CurrentPoint2D = FVector2D(PathPoints[i].X, PathPoints[i].Y);
		FVector2D NextPoint2D = FVector2D(PathPoints[(i + 1) % PathPoints.Num()].X, PathPoints[(i + 1) % PathPoints.Num()].Y);

		FVector CurrentPoint(CurrentPoint2D.X, CurrentPoint2D.Y, 0.0f);
		FVector NextPoint(NextPoint2D.X, NextPoint2D.Y, 0.0f);

		// Calculate the distance from the point to the segment
		float Distance = FMath::PointDistToSegment(Point, CurrentPoint, NextPoint);

		if (Distance <= (useOffset ? (Thickness + ThicknessOffset) : (Thickness + SmoothingThicknessOffset)))
		{
			return true;
		}
	}
	return false; // Point is not near the path
}


float AProceduralTerrain::DistFromPath(int32 X, int32 Y, bool useOffset) const
{
	FVector2D Point2D(X * Scale, Y * Scale);
	FVector Point(Point2D.X, Point2D.Y, 0.0f);
	float FinalDistance = INFINITY;
	for (int32 i = 0; i < PathPoints.Num(); ++i)
	{
		// Get current and next point, converting to FVector
		FVector2D CurrentPoint2D = FVector2D(PathPoints[i].X, PathPoints[i].Y);
		FVector2D NextPoint2D = FVector2D(PathPoints[(i + 1) % PathPoints.Num()].X, PathPoints[(i + 1) % PathPoints.Num()].Y);

		FVector CurrentPoint(CurrentPoint2D.X, CurrentPoint2D.Y, 0.0f);
		FVector NextPoint(NextPoint2D.X, NextPoint2D.Y, 0.0f);

		// Calculate the distance from the point to the segment
		float Distance = FMath::PointDistToSegment(Point, CurrentPoint, NextPoint);

		if (Distance <= FinalDistance)
		{
			FinalDistance = Distance;
		}
	}

	return (FinalDistance - (Thickness + ThicknessOffset));

}


FVector2D CatmullRomInterpolate(const FVector2D& P0, const FVector2D& P1, const FVector2D& P2, const FVector2D& P3, float T)
{
	float T2 = T * T;
	float T3 = T2 * T;

	return 0.5f * (
		(2.0f * P1) +
		(-P0 + P2) * T +
		(2.0f * P0 - 5.0f * P1 + 4.0f * P2 - P3) * T2 +
		(-P0 + 3.0f * P1 - 3.0f * P2 + P3) * T3
		);
}


void AProceduralTerrain::GeneratePath()
{
	// PathPoints.Empty();
	FRandomStream RandomStream(PathSeed);

	FVector2D CurrentPosition(-Width * Scale * 10.0f, Height * Scale * 0.5f); // Start position
	FVector2D CurrentDirection(1.0f, 0.0f); // Initial direction (X-axis)
	float StepSize = Scale * 5.0f; // Step distance per point

	TArray<FVector2D> BasePoints;
	BasePoints.Add(CurrentPosition);

	float MaxTurnAngle = 2.5f;     // Maximum allowed turn per section (in degrees)
	float CurrentTurnAngle = 0.0f;  // Keeps track of turn progress

	for (int32 i = 1; i < NumPoints; ++i)
	{
		// Compute a new target turn every `TurnSize` points
		if (i % TurnSize == 0)
		{
			// Generate a small, controlled turn angle (positive or negative)
			float NewTurnAngle = RandomStream.FRandRange(-MaxTurnAngle, MaxTurnAngle);
			CurrentTurnAngle = NewTurnAngle / TurnSmoothPoints; // Spread it over multiple points
		}

		// Apply the gradual turn
		float AngleRad = FMath::DegreesToRadians(CurrentTurnAngle);
		float CosA = FMath::Cos(AngleRad);
		float SinA = FMath::Sin(AngleRad);

		FVector2D NewDirection = FVector2D(
			CurrentDirection.X * CosA - CurrentDirection.Y * SinA,
			CurrentDirection.X * SinA + CurrentDirection.Y * CosA
		).GetSafeNormal();

		CurrentDirection = NewDirection;
		CurrentPosition += CurrentDirection * StepSize;
		BasePoints.Add(CurrentPosition);
	}

	// Apply Catmull-Rom interpolation to smooth the path further
	for (int32 i = 0; i < BasePoints.Num(); ++i)
	{
		const FVector2D& P0 = BasePoints[FMath::Max(0, i - 1)];
		const FVector2D& P1 = BasePoints[i];
		const FVector2D& P2 = BasePoints[FMath::Min(i + 1, BasePoints.Num() - 1)];
		const FVector2D& P3 = BasePoints[FMath::Min(i + 2, BasePoints.Num() - 1)];

		float T = 1.0f;
		FVector2D SmoothedPoint = CatmullRomInterpolate(P0, P1, P2, P3, T);
		FVector FinalPoint = FVector(SmoothedPoint.X, SmoothedPoint.Y, CalculateNoiseAtPoint(SmoothedPoint.X / Scale, SmoothedPoint.Y / Scale));
		if (NumPoints == 800) {
			PathPoints.Add(FinalPoint);
		}

	}
	NumPoints += 40;
}


void AProceduralTerrain::SmoothPathPointsHeight(float smoothLevel)
{
	if (NumPoints == 840) {
		for (int32 i = 1; i < PathPoints.Num(); ++i)
		{
			FVector& CurrentPoint = PathPoints[i];
			FVector& PrevPoint = PathPoints[i - 1];

			float InterpolatedHeight = FMath::Lerp(CurrentPoint.Z, PrevPoint.Z, smoothLevel);
			CurrentPoint.Z = InterpolatedHeight;
		}
	}
}


void AProceduralTerrain::GenerateTerrainSection(TerrainComponent* Component)
{
	// Launch async task for heavy computation
	Async(EAsyncExecution::ThreadPool, [=, this]()
		{
			////////////////////////////////////

			FRealtimeMeshStreamSet StreamSet;
			TRealtimeMeshBuilderLocal<uint16, FPackedNormal, FVector2DHalf, 1> Builder(StreamSet);

			// here we go ahead and enable all the basic mesh data parts
			Builder.EnableTangents();
			Builder.EnableTexCoords();
			Builder.EnablePolyGroups();


			int32 SectionSize = Width;
			int32 SectionX = Component->GetGridPosition().X;
			int32 SectionY = Component->GetGridPosition().Y;
			int32 StartX = SectionX * SectionSize;
			int32 StartY = SectionY * SectionSize;
			int32 EndX = StartX + SectionSize;
			int32 EndY = StartY + SectionSize;
			int32 LOD = pow(2, Component->GetLOD());
			int32 vertsPerRow = (SectionSize / LOD) + 1;
			int32 totalSize = vertsPerRow * vertsPerRow;



			// Generate vertices and UVs
			int vertsAmount = 0;

			for (int32 y = StartY; y <= EndY; y += LOD)
			{
				for (int32 x = StartX; x <= EndX; x += LOD)
				{
					float Z = CalculateHeight(x, y);
					Builder.AddVertex(FVector3f(x * Scale, y * Scale, Z))
						.SetTexCoord(FVector2f(
							static_cast<float>(x - StartX) / SectionSize,
							static_cast<float>(y - StartY) / SectionSize
						) * UVScale);
					vertsAmount++;
				}
			}

			int triangleCount = 0;
			for (int32 y = 0; y < SectionSize; y += LOD)
			{
				for (int32 x = 0; x < SectionSize; x += LOD)
				{
					int32 CurrentRow = y / LOD;
					int32 NextRow = (y + LOD) / LOD;

					int32 CurrentIndex = CurrentRow * vertsPerRow + (x / LOD);
					int32 RightIndex = CurrentIndex + 1;
					int32 BottomIndex = NextRow * vertsPerRow + (x / LOD);
					int32 BottomRightIndex = BottomIndex + 1;

					Builder.AddTriangle(CurrentIndex, BottomIndex, RightIndex, 0);
					Builder.AddTriangle(RightIndex, BottomIndex, BottomRightIndex, 0);
					triangleCount += 2;
				}
			}


			// Calculate normals and tangents without using texture coordinates
			TArray<FVector3f> VertexNormals;
			VertexNormals.SetNum(totalSize, false);
			TArray<FVector3f> VertexTangents;
			VertexTangents.SetNum(totalSize, false);
			for (int32 i = 0; i < totalSize; i++) {
				VertexNormals[i] = FVector3f::ZeroVector;
				VertexTangents[i] = FVector3f::ZeroVector;
			}

			for (int32 i = 0; i < triangleCount; i++)
			{
	

				TIndex3<uint32> Index = Builder.GetTriangle(i);
				const FVector3f& Vertex0 = Builder.GetPosition(Index[0]);
				const FVector3f& Vertex1 = Builder.GetPosition(Index[1]);
				const FVector3f& Vertex2 = Builder.GetPosition(Index[2]);

				// Calculate edges
				FVector3f Edge1 = Vertex1 - Vertex0;
				FVector3f Edge2 = Vertex2 - Vertex0;

				// Calculate the face normal using the cross product
				FVector3f FaceNormal = FVector3f::CrossProduct(Edge1, Edge2).GetSafeNormal() * -1.0f;

				// Accumulate normals for each vertex in the triangle
				VertexNormals[Index[0]] += FaceNormal;
				VertexNormals[Index[1]] += FaceNormal;
				VertexNormals[Index[2]] += FaceNormal;

				// Generate a tangent using an arbitrary orthogonal vector to the normal
				FVector3f Tangent = FVector3f::CrossProduct(FaceNormal, Edge1).GetSafeNormal();

				// Accumulate tangents for each vertex in the triangle
				VertexTangents[Index[0]] += Tangent;
				VertexTangents[Index[1]] += Tangent;
				VertexTangents[Index[2]] += Tangent;
			}

			// Normalize and set normals and tangents
			for (int32 i = 0; i < totalSize; i++)
			{
				VertexNormals[i].Normalize();
				//UE_LOG(LogTemp, Warning, TEXT("Normal %d: %s"), i, *VertexNormals[i].ToString());
				Builder.SetNormal(i, VertexNormals[i]);

				VertexTangents[i].Normalize();
				Builder.SetTangent(i, VertexTangents[i]);
			}


			Component->SetIsActive(true);
			if (!Component->GetIsInitialised()) {
				FString ComponentName = FString::Printf(TEXT("MainSection %i"), Component->GetIndex());
				FString ComponentName2 = FString::Printf(TEXT("MeshSection %i"), Component->GetIndex());
				FRealtimeMeshLODKey keyLOD = FRealtimeMeshLODKey::FRealtimeMeshLODKey(0);
				FRealtimeMeshSectionGroupKey GroupKey = FRealtimeMeshSectionGroupKey::Create(keyLOD, FName(ComponentName));
				RealtimeMesh->CreateSectionGroup(GroupKey, StreamSet);
				const FRealtimeMeshSectionKey Key = FRealtimeMeshSectionKey::Create(GroupKey, FName(ComponentName2));
				Component->SetGroupKey(GroupKey);
				Component->SetKey(Key);
				const FRealtimeMeshStreamRange StreamRange(0, totalSize - 1, 0, totalSize - 1);
				FRealtimeMeshSectionConfig sectionCongig;
				sectionCongig.bCastsShadow = false;
				sectionCongig.MaterialSlot = 0;
				sectionCongig.DrawType = ERealtimeMeshSectionDrawType::Static;
				bool hasCollision = Component->GetLOD() == 0;
				RealtimeMesh->CreateSection(Key, sectionCongig, StreamRange, hasCollision);
				RealtimeMesh->UpdateSectionConfig(Key, sectionCongig, hasCollision);
				RealtimeMesh->UpdateSectionGroup(GroupKey, StreamSet);
				Component->SetIsInitialised(true);
			}
			else if (RealtimeMesh != nullptr) {
				FRealtimeMeshSectionConfig sectionCongig;
				sectionCongig.bCastsShadow = false;
				sectionCongig.MaterialSlot = 0;
				sectionCongig.DrawType = ERealtimeMeshSectionDrawType::Static;
				bool hasCollision = Component->GetLOD() == 0;
				RealtimeMesh->UpdateSectionConfig(Component->GetKey(), sectionCongig, hasCollision);
				RealtimeMesh->UpdateSectionGroup(Component->GetGroupKey(), StreamSet);
			}
		});

}

void AProceduralTerrain::GeneratePathMesh()
{

	PathVertices.Empty();
	PathTriangles.Empty();
	PathNormals.Empty();
	PathUVs.Empty();

	FVector PreviousForward = FVector::ForwardVector;
	FVector LastLeftVertex;
	FVector LastRightVertex;
	TArray<FVector> LastVertexes;

	for (int32 i = 0; i < PathPoints.Num() - 1; ++i)
	{
		FVector StartPoint = PathPoints[i];
		FVector EndPoint = PathPoints[(i + 1)];
		float StartHeightCenter = StartPoint.Z;
		float EndHeightCenter = EndPoint.Z;
		FVector CurrentPoint3D = FVector(StartPoint.X, StartPoint.Y, StartHeightCenter + HeightAdjust);
		FVector NextPoint3D = FVector(EndPoint.X, EndPoint.Y, EndHeightCenter + HeightAdjust);
		FVector Forward = (NextPoint3D - CurrentPoint3D).GetSafeNormal();
		FVector Up = FVector::UpVector;
		FVector Right = FVector::CrossProduct(Up, Forward).GetSafeNormal();

		TArray<FVector> CurrentNextVertexes;
		int32 IndexOffset = PathVertices.Num();
		for (int l = 0; l < (ThicknessDetail * 2) + 1; l++) {
			if (l == ThicknessDetail) {
				if (i == 0) {
					PathVertices.Add(NextPoint3D);
					PathVertices.Add(CurrentPoint3D);
					CurrentNextVertexes.Add(NextPoint3D);
				}
				else {
					PathVertices.Add(NextPoint3D);
					PathVertices.Add(LastVertexes[l]);
					CurrentNextVertexes.Add(NextPoint3D);
				}
			}
			else if (l < ThicknessDetail) {
				// left
				float EndHeightLeft = CalculateNoiseAtPoint((EndPoint.X - ((Thickness / ThicknessDetail) * (ThicknessDetail - l)) * Right.X) / Scale, (EndPoint.Y - ((Thickness / ThicknessDetail) * (ThicknessDetail - l)) * Right.Y) / Scale);
				EndHeightLeft = FMath::Lerp(EndHeightLeft, EndHeightCenter, Flatness);
				FVector EndLeft = FVector(EndPoint.X - ((Thickness / ThicknessDetail) * (ThicknessDetail - l)) * Right.X, EndPoint.Y - ((Thickness / ThicknessDetail) * (ThicknessDetail - l)) * Right.Y, EndHeightLeft + HeightAdjust);
				CurrentNextVertexes.Add(EndLeft);
				if (i == 0) {
					float StartHeightLeft = CalculateNoiseAtPoint((StartPoint.X - ((Thickness / ThicknessDetail) * (ThicknessDetail - l)) * Right.X) / Scale, (StartPoint.Y - ((Thickness / ThicknessDetail) * (ThicknessDetail - l)) * Right.Y) / Scale);
					StartHeightLeft = FMath::Lerp(StartHeightLeft, StartHeightCenter, Flatness);
					FVector StartLeft = FVector(StartPoint.X - ((Thickness / ThicknessDetail) * (ThicknessDetail - l)) * Right.X, StartPoint.Y - ((Thickness / ThicknessDetail) * (ThicknessDetail - l)) * Right.Y, StartHeightLeft + HeightAdjust);

					PathVertices.Add(EndLeft);
					PathVertices.Add(StartLeft);
				}
				else {
					PathVertices.Add(EndLeft);
					PathVertices.Add(LastVertexes[l]);
				}
			}
			else {
				//right
				float EndHeightRight = CalculateNoiseAtPoint((EndPoint.X + ((Thickness / ThicknessDetail) * (l - ThicknessDetail)) * Right.X) / Scale, (EndPoint.Y + ((Thickness / ThicknessDetail) * (l - ThicknessDetail)) * Right.Y) / Scale);
				EndHeightRight = FMath::Lerp(EndHeightRight, EndHeightCenter, Flatness);
				FVector EndRight = FVector(EndPoint.X + ((Thickness / ThicknessDetail) * (l - ThicknessDetail)) * Right.X, EndPoint.Y + ((Thickness / ThicknessDetail) * (l - ThicknessDetail)) * Right.Y, EndHeightRight + HeightAdjust);
				CurrentNextVertexes.Add(EndRight);
				if (i == 0) {

					float StartHeightRight = CalculateNoiseAtPoint((StartPoint.X + ((Thickness / ThicknessDetail) * (l - ThicknessDetail)) * Right.X) / Scale, (StartPoint.Y + ((Thickness / ThicknessDetail) * (l - ThicknessDetail)) * Right.Y) / Scale);
					StartHeightRight = FMath::Lerp(StartHeightRight, StartHeightCenter, Flatness);
					FVector StartRight = FVector(StartPoint.X + ((Thickness / ThicknessDetail) * (l - ThicknessDetail)) * Right.X, StartPoint.Y + ((Thickness / ThicknessDetail) * (l - ThicknessDetail)) * Right.Y, StartHeightRight + HeightAdjust);

					PathVertices.Add(EndRight);
					PathVertices.Add(StartRight);
				}
				else {
					PathVertices.Add(EndRight);
					PathVertices.Add(LastVertexes[l]);
				}
			}
		}
		LastVertexes = CurrentNextVertexes;

		for (int j = 0; j < (ThicknessDetail * 2); j++) {
			PathTriangles.Add(IndexOffset + 1 + (j * 2));
			PathTriangles.Add(IndexOffset + 2 + (j * 2));
			PathTriangles.Add(IndexOffset + 0 + (j * 2));

			PathTriangles.Add(IndexOffset + 3 + (j * 2));
			PathTriangles.Add(IndexOffset + 2 + (j * 2));
			PathTriangles.Add(IndexOffset + 1 + (j * 2));
		}

		// UV coordinates
		for (int l = 0; l < (ThicknessDetail * 2) + 1; l++) {
			float UCoord = static_cast<float>(l) / (ThicknessDetail * 2);
			float VCoord = i * PathTextureScale;
			PathUVs.Add(FVector2D(UCoord, VCoord));
			PathUVs.Add(FVector2D(UCoord, (i + 1) * PathTextureScale));
		}
	}

}



void AProceduralTerrain::DisplayPathMesh() {

	for (int i = 0; i < PathVertices.Num(); i) {
		PathVertices[i].Z += EdgeHeightOffset;
		PathVertices[i + 1].Z += EdgeHeightOffset;

		PathVertices[i + (((ThicknessDetail * 2) + 1) * 2) - 2].Z += EdgeHeightOffset;
		PathVertices[i + (((ThicknessDetail * 2) + 1) * 2) - 1].Z += EdgeHeightOffset;

		i += ((ThicknessDetail * 2) + 1) * 2;
	}


	// Mesh
	FString ComponentName = FString::Printf(TEXT("Path"));
	UProceduralMeshComponent* CurrentMeshComponent = NewObject<UProceduralMeshComponent>(this, FName(*ComponentName), RF_Transactional);
	MeshComponent.Add(CurrentMeshComponent);
	CurrentMeshComponent->SetupAttachment(RootComponent);
	CurrentMeshComponent->RegisterComponent();
	CurrentMeshComponent->CreateMeshSection_LinearColor(0, PathVertices, PathTriangles, PathNormals, PathUVs, TArray<FLinearColor>(), TArray<FProcMeshTangent>(), true);
	CurrentMeshComponent->SetMaterial(0, PathMaterial);
	CurrentMeshComponent->ContainsPhysicsTriMeshData(true);

}
